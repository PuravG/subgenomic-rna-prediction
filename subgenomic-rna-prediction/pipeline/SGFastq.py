#!/usr/bin/env python3
import os
import gzip
import shutil
import pathlib
import SGErrors


def zip_to_file(infile, outfile):
    with open(infile, 'rb') as f_in:
        with gzip.open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out, )


def zip_and_link(fastq_s, fastq_p1, fastq_p2, output):
    outfile_s, outfile_p1, outfile_p2 = output
    if fastq_s:
        if os.path.exists(fastq_s):
            if fastq_s.endswith(".gz"):
                ap = os.path.abspath(fastq_s)
                os.symlink(ap, outfile_s)

            else:
                zip_to_file(fastq_s, outfile_s)
            pathlib.Path(outfile_p1).touch(exist_ok=True)
            pathlib.Path(outfile_p2).touch(exist_ok=True)
        else:
            SGErrors.fastq_not_found([fastq_s])

    elif fastq_p1 and fastq_p2:
        if os.path.exists(fastq_p1) and os.path.exists(fastq_p2):
            if fastq_p1.endswith(".gz"):
                ap1 = os.path.abspath(fastq_p1)
                ap2 = os.path.abspath(fastq_p2)
                os.symlink(ap1, outfile_p1)
                os.symlink(ap2, outfile_p2)
            else:
                zip_to_file(fastq_p1, outfile_p1)
                zip_to_file(fastq_p2, outfile_p2)
            pathlib.Path(outfile_s).touch(exist_ok=True)
        else:
            SGErrors.fastq_not_found([fastq_p1, fastq_p2])
    else:
        SGErrors.no_files_found([fastq_s, fastq_p1, fastq_p2])


def touch_only(outfiles):
    pass
