from Bio import SeqIO

def extract_rDNA_read_ids(gff_files):
    read_ids = set()
    for gff_file in gff_files:
        with open(gff_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")  # explicit tab split for GFF
                read_ids.add(fields[0].strip())  # strip any stray whitespace
    return read_ids


def filter_fasta(input_fasta, output_fasta, read_ids):
    count = 0
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id.strip() in read_ids:
                SeqIO.write(record, out_f, "fasta")
                count += 1
    return count


def run(ctx):
    filtered = ctx.artifact("rdna_reads", "reads", ".fasta")
    
    if ctx.artifact_exists_or_skip("rdna_reads"):
        return

    # fasta is the converted FASTA. input_fasta is the original file, which may be FASTQ.
    input_fasta = ctx.require("fasta")
    bac = ctx.require("barrnap_bac")
    euk = ctx.require("barrnap_euk")
 
    read_ids = extract_rDNA_read_ids([bac, euk])
    ctx.log(f"Input FASTA exists: {input_fasta.exists()}, path: {input_fasta}")
    ctx.log(f"BAC GFF exists: {bac.exists()}, path: {bac}")
    ctx.log(f"EUK GFF exists: {euk.exists()}, path: {euk}")
    ctx.log(f"Number of read IDs extracted: {len(read_ids)}")
    ctx.log(f"Some IDs: {list(read_ids)[:5]}")

    read_ids = extract_rDNA_read_ids([bac, euk])
    if not read_ids:
        ctx.log("No rDNA reads found — terminating pipeline early.")
        ctx.register("pipeline_terminated_early", True)
        return

    ctx.log(f"Found {len(read_ids)} rDNA hits; filtering input FASTA...")
    written = filter_fasta(input_fasta, filtered, read_ids)
    ctx.log(f"Wrote {written} sequences to {filtered}")