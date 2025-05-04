import subprocess as sp
import logging

def run_cmd(cmd: list[str], shell: bool=False, silence: str=False):
    logger = logging.getLogger("TDNAreader")
    logger.info(f"Running: {cmd}\n")
    if silence:
        sp.run(cmd, shell=True, stdout=silence, stderr=silence)
    else:
        sp.run(cmd, shell=True)