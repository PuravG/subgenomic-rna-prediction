#!/usr/bin/env python3

def fastq_not_found(filenames):
    raise RuntimeError("""
###
At least one of your FASTQ files does not exist: %s
###
""" % (",".join(filenames)))


def no_files_found(filenames):
    raise RuntimeError("""
###
FASTQ input specified but no fastq files provided
###
""")


def bad_start_point(start_point):
    raise RuntimeError(
                """
###
The start_point parameter should be sra_id, raw_fastq, trimmed_fastq or bam, \
not %s
###
""" % start_point)
