#!/usr/bin/env python3
import logging
import SGErrors
import subprocess


def log_command(cmdout, command):
    cmd = open(cmdout, "a")
    cmd.write(f"{command}\n")
    cmd.close()


def set_up_logging(outdir, prefix):
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


def write_check_log(captured: str,
                    outdir: str,
                    stem: str,
                    ID: str,
                    normal_length: int,
                    type: str):
    # Output log info
    if type == 'stdout':
        log = captured.stdout.decode().split("\n")
    elif type == 'stderr':
        log = captured.stderr.decode().split("\n")
    out = open(f"{outdir}/{stem}.log", "w")
    out.write("\n".join(log))
    out.close()
    # If it worked the log should be 4 lines
    if len(log) != normal_length:
        SGErrors.tool_fail(stem, ID)


def touch_only(outfiles, cmdout, logi):
    for outfile in outfiles:
        statement = f'touch {outfile}'
        log_command(cmdout, statement)
        subprocess.run(statement, shell=True)
