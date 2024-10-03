#!/usr/bin/env python3
import os
import pathlib
import SGErrors
import subprocess


def zip_to_file(infile, outfile):
    '''
    Compress file using pigz and save the output to outfile

    Parameters
    ----------
    infile: str
        Path to input file
    outfile: str
        Path to desired output file

    Returns
    -------
    None
    '''
    # Links need an absolute path
    ap = os.path.abspath(infile)
    statement = f'pigz -c {ap} > {outfile}'
    subprocess.run(statement, shell=True)


def zip_and_link(fastq_s, fastq_p1, fastq_p2, outfiles):
    '''
    fastq_s is the path to an input single end FASTQ file, fastq_p1 and
    fastq_p2 to input paired end fastq files - only either fastq_s or 
    fastq_p1 and fastq_p2 need to be valid file paths.

    If the fastq files are already gzipped (assumed if suffix is .gz), just
    link to them.
    If they are not zipped, zip with pigz.
    Make empty placeholder files for single end if the input is paired end
    and vice versa.

    Parameters
    ----------
    fastq_s: str
        Path to single end fastq file, if this is a valid file path it will
        be used and the paired end files ignored
    fastq_p1:
        Path to read one fastq file for paired end
    fastq_p2:
        Path to read two fastq file for paired end
    outfiles: list
        List of paths to single end, paired end 1, paired end 2 desired output
        files

    Returns
    -------
    None
    '''
    outfile_s, outfile_p1, outfile_p2 = outfiles
    # Always use the single end file if there is one
    if fastq_s:
        # Single end
        if os.path.exists(fastq_s):
            # Assume ending with gz = gzipped
            if fastq_s.endswith(".gz"):
                # link to the file
                ap = os.path.abspath(fastq_s)
                os.symlink(ap, outfile_s)
            else:
                # gzip the file
                zip_to_file(fastq_s, outfile_s)
            # Empty placeholders for paired end
            pathlib.Path(outfile_p1).touch(exist_ok=True)
            pathlib.Path(outfile_p2).touch(exist_ok=True)
        else:
            # If the file path is not valid
            SGErrors.fastq_not_found([fastq_s])
    elif fastq_p1 and fastq_p2:
        # Paired end
        if os.path.exists(fastq_p1) and os.path.exists(fastq_p2):
            # Assume ending with gz = gzipped
            if fastq_p1.endswith(".gz"):
                # Links need the absolute path
                ap1 = os.path.abspath(fastq_p1)
                ap2 = os.path.abspath(fastq_p2)
                # Link to the files
                os.symlink(ap1, outfile_p1)
                os.symlink(ap2, outfile_p2)
            else:
                # Gzip the files
                zip_to_file(fastq_p1, outfile_p1)
                zip_to_file(fastq_p2, outfile_p2)
            # Empty placeholder for single end
            pathlib.Path(outfile_s).touch(exist_ok=True)
        else:
            # If the file path is not valid
            SGErrors.fastq_not_found([fastq_p1, fastq_p2])
    else:
        # If both file paths are None or False raise an error
        SGErrors.no_files_found([fastq_s, fastq_p1, fastq_p2])


def download_fastq(sra_id, method, outfiles, outdir):
    '''
    Download a FASTQ file using the specified method.

    Parameters
    ----------
    sra_id: str
        SRA run ID to download
    method: str
        Can be sra, ena or ena_aspera to download with NCBI fasterq-dump,
        enaDataGet or enaDataGet plus Aspera, respectively.
    outfiles: list
        List of paths to single end, paired end 1, paired end 2 desired output
        files
    outdir: str
        Output directory for log files etc

    Returns
    -------
    None
    '''
    if method == 'sra':
        statement = f'fasterq-dump --split-files {sra_id} --outdir {outdir}'
        fasterq = subprocess.run(statement, shell=True, capture_output=True)
        log = fasterq.stderr.decode().split("\n")
        out = open(f"{outdir}/fasterq_dump.log", "w")
        out.write("\n".join(log))
        out.close()
        if len(log) != 4:
            SGErrors.fasterq_fail(sra_id)
        if os.path.exists(f'{outdir}/{sra_id}_1.fastq') and os.path.exists(
                f'{outdir}/{sra_id}_2.fastq'):
            zip_to_file(f'{outdir}/{sra_id}_1.fastq', outfiles[1])
            zip_to_file(f'{outdir}/{sra_id}_2.fastq', outfiles[2])
            pathlib.Path(outfiles[0]).touch(exist_ok=True)
            os.remove(f'{outdir}/{sra_id}_1.fastq')
            os.remove(f'{outdir}/{sra_id}_2.fastq')
        elif os.path.exists(f'{outdir}/{sra_id}.fastq'):
            zip_to_file(f'{outdir}/{sra_id}.fastq', outfiles[0])
            pathlib.Path(outfiles[1]).touch(exist_ok=True)
            pathlib.Path(outfiles[1]).touch(exist_ok=True)
            os.remove(f'{outdir}/{sra_id}.fastq')


def touch_only(outfiles):
    for outfile in outfiles:
        pathlib.Path(outfile).touch(exist_ok=True)
