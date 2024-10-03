#!/usr/bin/env python3
import logging
import SGErrors
import subprocess
from typing import Optional
from typing import TextIO


def set_up_logging(outdir: str,
                   prefix: str):
    '''
    Set up a logger opject which writes to both STDOUT and to file.

    Parameters
    ----------
    outdir: str
        Directory to store the log file
    prefix: str
        Name for the log file - it will be saved as outdir/prefix.log

    Returns
    -------
    logi: logging.Logger
        An open log file object
    '''
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


def run_log_command(command: str,
                    cmdout: str,
                    outdir: Optional[str] = None,
                    toolname: Optional[str] = None,
                    ID: Optional[str] = None,
                    logi: Optional[TextIO] = None):
    '''
    For a bash command - write the command to the bash command log,
    run the command and store the output, if the other parameters are
    specified then write the stdout and stdlog from the command to file
    and check them for errors, check the return code was 0.
    Parameters
    ----------
    command: str
        A bash statement to pass to the command line
    cmdout: str
        Path to file to output the bash command
    outdir: str
        If specified, copy stdout and stderr to this directory
    toolname: str
        Name of the tool - to give informative filenames
    ID: str
        Sample ID - for informative filenames
    logi: logging.Logger
        Open log file
    Returns
    -------
    cap:
        The captured output of the bash command
    '''
    # Store the command
    cmd = open(cmdout, "a")
    cmd.write(f"{command}\n")
    cmd.close()
    if logi:
        logi.info(f"Running command: {command}")
    # Run the command
    cap = subprocess.run(command, shell=True, capture_output=True)

    if outdir:
        # Output STDOUT and STDERR to file
        sout = cap.stdout.decode().split("\n")
        serr = cap.stderr.decode().split("\n")

        out = open(f"{outdir}/{toolname}_{ID}_stderr.log", "w")
        out.write("\n".join(serr))
        out.close()
        out = open(f"{outdir}/{toolname}_{ID}_stdout.log", "w")
        out.write("\n".join(sout))
        out.close()
    if cap.returncode != 0:
        SGErrors.tool_fail(toolname, ID)

    return (cap)


def touch_only(outfiles: list,
               cmdout: str,
               logi: TextIO):
    '''
    Uses the Linux touch command to make empty placeholder output files
    where they do not already exist.

    Parameters
    ----------
    outfiles: list
        List of files to touch
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file

    Returns
    -------
    None
    '''
    for outfile in outfiles:
        statement = f'touch {outfile}'
        run_log_command(statement, cmdout)
