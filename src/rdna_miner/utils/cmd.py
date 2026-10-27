import subprocess

from rdna_miner.utils.logging_utils import info

def run_command(command: str, log_file=None):
    info(f"Running: {command}")

    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True
    )

    if log_file:
        with open(log_file, "a") as log:
            log.write(result.stdout)
            log.write(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{command}\n{result.stderr}")

    return result