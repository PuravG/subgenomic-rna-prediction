#!/usr/bin/env python3

def fastq_not_found(filenames):
    fn = ", ".join(filenames)
    raise RuntimeError(f"""
###
At least one of your FASTQ files does not exist: {fn}
###
""")


def no_files_found():
    raise RuntimeError("""
###
FASTQ input specified but no fastq files provided
###
""")


def bad_start_point(start_point):
    raise RuntimeError(f"""
###
The start_point parameter should be sra_id, raw_fastq, trimmed_fastq or bam, \
not {start_point}
###
""")


def tool_fail(stem, sra_id):
    raise RuntimeError(
                f"""
###
{stem} run failed with ID {sra_id}
###
""")


def bad_download_source(method):
    raise RuntimeError(
                f"""
###
Download method {method} not found
###
""")


def no_aspera_config(aspera_config):
    raise FileNotFoundError(
                f"""
###
Aspera config file not found at {aspera_config}. Specify a valid file or
use source "ena" instead of "ena_aspera" to download without Aspera
###
""")


def trimmer_not_found(trimmer):
    raise RuntimeError(
                f"""
###
Trimming tool {trimmer} is not currently implemented
###
""")
