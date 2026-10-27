from rdna_miner.steps import (barrnap_scan,
                              filter_rdna_reads,
                              assemble_rdna,
                              annotate_rfam,
                              taxonomy_decipher,
                              map_reads,
                              compile_taxonomy,)
from rdna_miner.utils.logging_utils import section


def build_long_read_pipeline():
    """
    Each step exposes run(ctx) and internally
    handles its own dependencies using ctx.require().
    """

    return [
        (barrnap_scan.run, "Running Barrnap"),
        (filter_rdna_reads.run, "Filtering putative rDNA reads"),
        (assemble_rdna.run, "Assembling rDNA reads with Flye"),
        (annotate_rfam.run, "Annotating rDNA operons with Rfam/cmscan"),
        (taxonomy_decipher.run, "Assigning taxonomy with DECIPHER"),
        (map_reads.run, "Mapping reads to rDNA assembly"),
        (compile_taxonomy.run, "Completing read-level taxonomy assignment"),
    ]

def build_assembly_pipeline():
    """Pipeline provided an assembly."""

    return [
        (annotate_rfam.run, "Annotating rDNA operons with Rfam/cmscan"),
        (taxonomy_decipher.run, "Assigning taxonomy with DECIPHER"),
    ]

def run_pipeline(ctx, pipeline):
    total_steps = len(pipeline)

    # The outer timer is the parent of each step row. Child-process RSS
    # inside a step covers Barrnap, Flye, cmscan, and R.
    with ctx.benchmark_timer("pipeline", input_file=ctx.input_file):
        for i, (step_func, step_name) in enumerate(pipeline, start=1):
            if ctx.exists("pipeline_terminated_early"):
                return

            section(f"[{i}/{total_steps}] {step_name}")
            with ctx.benchmark_timer(step_name, input_file=ctx.input_file):
                step_func(ctx)

        ctx.report_artifacts()