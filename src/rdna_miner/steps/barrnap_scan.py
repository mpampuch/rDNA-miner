import re
from pathlib import Path
from rdna_miner.utils.cmd import run_command

# 1.10.x dropped --kingdom euk and --reject. This pipeline still calls both.
_BARRNAP_VERSION = re.compile(r"barrnap\s+(\d+\.\d+(?:\.\d+)?)", re.IGNORECASE)
_REQUIRED_BARRNAP_VERSION = "0.9"


def barrnap_version_is_supported(version_text: str) -> bool:
    """True only for Barrnap 0.9, the last release with --kingdom euk and --reject."""
    match = _BARRNAP_VERSION.search(version_text)
    return match is not None and match.group(1) == _REQUIRED_BARRNAP_VERSION


def require_barrnap_0_9(version_text: str) -> None:
    if barrnap_version_is_supported(version_text):
        return
    reported = version_text.strip() or "unknown"
    raise RuntimeError(
        "Barrnap 0.9 is required (--kingdom euk and --reject were removed in 1.10.x). "
        f"Reported version was: {reported}"
    )

def gff_has_records(gff_file: Path) -> bool:
    """Return True if GFF has any non-comment lines besides the header."""
    with open(gff_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            return True
    return False

def run(ctx):
    bac_gff = ctx.artifact("barrnap_bac", "barrnap", ".gff")
    euk_gff = ctx.artifact("barrnap_euk", "barrnap", ".gff")
    combined_gff = ctx.artifact("barrnap_combined", "barrnap", ".gff")

    if ctx.artifact_exists_or_skip("barrnap_combined"):
        return

    version = run_command("barrnap --version")
    require_barrnap_0_9(f"{version.stdout}\n{version.stderr}")

    ctx.log(f"Running Barrnap for bacteria: {ctx.fasta} -> {bac_gff}")
    run_command(f"barrnap --kingdom bac --threads {ctx.threads} {ctx.fasta} --reject 0.2 > {bac_gff}")

    ctx.log(f"Running Barrnap for eukaryotes: {ctx.fasta} -> {euk_gff}")
    run_command(f"barrnap --kingdom euk --threads {ctx.threads} {ctx.fasta} --reject 0.2 > {euk_gff}")

    ctx.log(f"Combining GFFs -> {combined_gff}")
    run_command(f"cat {bac_gff} {euk_gff} > {combined_gff}")

    if not gff_has_records(combined_gff):
        ctx.terminate_pipeline("No rDNA hits found in the input FASTA")
        return