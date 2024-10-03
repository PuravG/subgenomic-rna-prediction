#!/usr/bin/env python3
import logging
import SGErrors
import subprocess
from typing import Optional
from typing import TextIO


def run_log_command(cmdout: str,
                    command: str,
                    outdir: Optional[str] = None,
                    toolname: Optional[str] = None,
                    ID: Optional[str] = None,
                    logi: Optional[TextIO] = None):
    '''
    Write a bash command to the command log
    '''
    cmd = open(cmdout, "a")
    cmd.write(f"{command}\n")
    cmd.close()
    cap = subprocess.run(command, shell=True, capture_output=True)
    x = 0
    if outdir:
        # Output log info
        sout = cap.stdout.decode().split("\n")
        serr = cap.stderr.decode().split("\n")

        out = open(f"{outdir}/{toolname}_{ID}_stderr.log", "w")
        out.write("\n".join(serr))
        out.close()
        out = open(f"{outdir}/{toolname}_{ID}_stdout.log", "w")
        out.write("\n".join(sout))
        out.close()

        for line in sout:
            if 'error' in line.lower():
                x += 1
        for line in serr:
            if 'error' in line.lower():
                x += 1
    # If it worked the log should be 4 lines
    if x != 0 or cap.returncode != 0:
        SGErrors.tool_fail(toolname, ID)

    return (cap)


def set_up_logging(outdir: str,
                   prefix: str):
    # Create handlers
    file_handler = logging.FileHandler(f'{outdir}/{prefix}.log')
    console_handler = logging.StreamHandler()

    # Create formatter and add it to handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Set up the root logger
    logi = logging.getLogger()
    logi.setLevel(logging.INFO)
    logi.addHandler(file_handler)
    logi.addHandler(console_handler)
    return (logi)


def touch_only(outfiles: list,
               cmdout: str,
               logi: TextIO):
    for outfile in outfiles:
        statement = f'touch {outfile}'
        run_log_command(cmdout, statement)
