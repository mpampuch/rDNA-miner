from pathlib import Path
from typing import Dict, Union, List
import logging
from rdna_miner.utils.logging_utils import section, info, console

from rdna_miner.utils.fasta import load_fasta_or_convert

ArtifactType = Union[Path, List[Path]]

class Context:
    """Holds input/output paths, artifacts, and the per-run benchmark log.

    script_name is written on every benchmark row so long-read and assembly
    runs can be told apart in the same TSV.
    """

    script_name = None

    SUBDIRS = {
        "summary": "summary",
        "reads": "reads",
        "assembly": "assembly",
        "taxonomy": "taxonomy",
        "mapping": "mapping",
        "figures": "figures"
    }

    def __init__(
        self,
        input_file: Path,
        output_dir: Path,
        force: bool = False,
        script_name: str | None = None,
    ):
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.force = force
        self.script_name = script_name
        self.artifacts: Dict[str, object] = {}  # Now can store any object

        self.logger = logging.getLogger("rdna-miner")

        # Timed because converting a large FASTQ can take minutes and would
        # otherwise be an unexplained gap before the first pipeline step.
        with self.benchmark_timer("fasta_convert", input_file=self.input_file):
            self.fasta = load_fasta_or_convert(self.input_file, self.output_dir)
        self.register("fasta", self.fasta)

        # Register original input fasta
        self.register("input_fasta", self.input_file)

    @property
    def benchmark_dir(self) -> Path:
        """Per-run benchmark log location: {output_dir}/benchmarking/."""
        d = self.output_dir / "benchmarking"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def benchmark_timer(self, label, **kwargs):
        from rdna_miner.utils.benchmark import BenchmarkTimer
        kwargs.setdefault("log_dir", str(self.benchmark_dir))
        if self.script_name and "script_name" not in kwargs:
            kwargs["script_name"] = self.script_name
        return BenchmarkTimer(label=label, **kwargs)

    # ---------------- Artifact Registration ----------------

    def register(self, name: str, value):
        """Register an artifact or runtime value. Accepts paths, lists, bools, or any object."""
        if isinstance(value, (list, tuple)):
            self.artifacts[name] = [Path(p) if not isinstance(p, Path) else p for p in value]
        elif isinstance(value, Path):
            self.artifacts[name] = value
        else:
            self.artifacts[name] = value  # can be bool, str, or any object

    def get(self, name: str):
        if name not in self.artifacts:
            raise ValueError(f"Artifact '{name}' not registered.")
        return self.artifacts[name]

    def require(self, name: str):
        artifact = self.get(name)

        if isinstance(artifact, list):
            for p in artifact:
                if not isinstance(p, Path) or not p.exists():
                    raise FileNotFoundError(p)
        elif isinstance(artifact, Path):
            if not artifact.exists():
                raise FileNotFoundError(artifact)

        return artifact

    # ---------------- Artifact Paths ----------------

    def artifact_path(self, name: str, subdir: str = None, suffix: str = None) -> Path:
        base = self.fasta.stem
        suffix = suffix or ".dat"

        path = self.output_dir
        if subdir:
            path = path / subdir
            path.mkdir(parents=True, exist_ok=True)

        return path / f"{base}_{name}{suffix}"

    def artifact(self, name: str, subdir: str = None, suffix: str = None) -> Path:
        path = self.artifact_path(name, subdir, suffix)
        self.register(name, path)
        return path

    # ---------------- Existence Checks ----------------

    def exists(self, name: str) -> bool:
        if self.force:
            return False

        artifact = self.artifacts.get(name)
        if artifact is None:
            return False

        if isinstance(artifact, list):
            return all(isinstance(p, Path) and p.exists() for p in artifact)
        elif isinstance(artifact, Path):
            return artifact.exists()

        # For non-path artifacts (like flags), existence is True if registered
        return True

    def artifact_exists_or_skip(self, name: str) -> bool:
        if name not in self.artifacts:
            return False

        exists = self.exists(name)
        if exists:
            info(f"Skipping '{name}' (artifact exists: {self.artifacts[name]})")
        return exists

    # ---------------- Logging ----------------

    def log_step(self, name: str):
        section(name)

    def log(self, message: str):
        info(message)

    # ---------------- Pipeline Helpers ----------------

    def terminate_pipeline(self, reason: str):
        """Set a flag to terminate the pipeline early and log reason."""
        self.log(f"Pipeline terminated early: {reason}")
        self.register("pipeline_terminated_early", True)

    def report_artifacts(self):
        """Print a nicely formatted summary of all registered artifacts."""
        console.print()
        console.rule("Generated output", align="left", style="Dim")
        for name, path in self.artifacts.items():
            if isinstance(path, list):
                info(f"{name}:")
                for p in path:
                    info(f"- {p}")
            else:
                info(f"{name}: {path}")
