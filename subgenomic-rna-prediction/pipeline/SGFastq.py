#!/usr/bin/env python3
import os
import SGErrors
import SGSetup
from typing import Optional
from typing import TextIO
import shutil


def zip_to_file(infile: str, outfile: str, cmdout: str, logi: TextIO):
    '''
    Compress file using pigz and save the output to outfile

    Parameters
    ----------
    infile: str
        Path to input file
    outfile: str
        Path to desired output file
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file

    Returns
    -------
    None
    '''
    logi.info(f"Zipping file {infile}")
    # Links need an absolute path
    ap = os.path.abspath(infile)
    statement = f'pigz -c {ap} > {outfile}'
    SGSetup.run_log_command(statement, cmdout)


def zip_and_link(fastq_s: Optional[str],
                 fastq_p1: Optional[str],
                 fastq_p2: Optional[str],
                 outfiles: list,
                 outdir: str,
                 cmdout: str, logi: TextIO):
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
    fastq_p1: str
        Path to read one fastq file for paired end
    fastq_p2: str
        Path to read two fastq file for paired end
    outfiles: list
        List of paths to single end, paired end 1, paired end 2 desired output
        files
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file
    Returns
    -------
    None
    '''

    outfile_s, outfile_p1, outfile_p2 = outfiles
    # Always use the single end file if there is one
    if fastq_s:
        logi.info(f"Single end FASTQ files provided: {fastq_s}")
        # Single end
        if os.path.exists(fastq_s):
            # Assume ending with gz = gzipped
            if fastq_s.endswith(".gz"):
                logi.info(f"FASTQ {fastq_s} is already zipped - creating link")
                # link to the file
                ap = os.path.abspath(fastq_s)
                statement = f'ln -s {ap} {outfile_s}'
                SGSetup.run_log_command(statement, cmdout)
            else:
                logi.info(f"FASTQ {fastq_s} is unzipped - zipping")
                # gzip the file
                zip_to_file(fastq_s, outfile_s, cmdout, logi)
            # Empty placeholders for paired end
            SGSetup.touch_only([outfile_p1, outfile_p2], cmdout, logi)
            typ = 'single'
        else:
            # If the file path is not valid
            SGErrors.fastq_not_found([fastq_s])
    elif fastq_p1 and fastq_p2:
        # Paired end
        if os.path.exists(fastq_p1) and os.path.exists(fastq_p2):
            # Assume ending with gz = gzipped
            if fastq_p1.endswith(".gz"):
                logi.info(
                    f"FASTQs {fastq_p1} and {fastq_p2} are already zipped - \
                      creating links")
                # Links need the absolute path
                ap1 = os.path.abspath(fastq_p1)
                ap2 = os.path.abspath(fastq_p2)
                # Link to the files
                os.symlink(ap1, outfile_p1)
                os.symlink(ap2, outfile_p2)
            else:
                # Gzip the files
                logi.info(f"FASTQ {fastq_p1} is unzipped - zipping")
                zip_to_file(fastq_p1, outfile_p1, cmdout, logi)
                logi.info(f"FASTQ {fastq_p2} is unzipped - zipping")
                zip_to_file(fastq_p2, outfile_p2, cmdout, logi)
            # Empty placeholder for single end
            SGSetup.touch_only([outfile_s], cmdout, logi)
            typ = 'paired'
        else:
            # If the file path is not valid
            SGErrors.fastq_not_found([fastq_p1, fastq_p2])
    else:
        # If both file paths are None or False raise an error
        SGErrors.no_files_found()
    typ_out = open(f"{outdir}/data_type.txt", "w")
    typ_out.write(typ)
    typ_out.close()


def download_fastq_fasterqdump(sra_id: str, outdir: str, outfiles: list,
                               cmdout: str, logi: TextIO):
    """
    Wrapper to download files using NCBI fasterq-dump.

    Parameters
    ----------
    sra_id: str
        SRA run ID to download
    outdir: str
        Output directory for log files etc
    outfiles: list
        List of paths to single end, paired end 1, paired end 2 desired output
        files
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file
    """
    # Download the files
    statement = f'fasterq-dump --split-files {sra_id} --outdir {outdir}'
    SGSetup.run_log_command(statement, cmdout, outdir, 'fasterq',
                            sra_id, logi)

    p1 = f'{outdir}/{sra_id}_1.fastq'
    p2 = f'{outdir}/{sra_id}_2.fastq'
    s = f'{outdir}/{sra_id}.fastq'
    # Determine paired end or single end, zip, remove unzipped
    if os.path.exists(p1) and os.path.exists(p2):
        zip_to_file(p1, outfiles[1], cmdout, logi)
        zip_to_file(p2, outfiles[2], cmdout, logi)
        os.remove(p1)
        os.remove(p2)
        return ("paired")
    else:
        zip_to_file(s, outfiles[0], cmdout, logi)
        os.remove(s)
        return ("single")


def download_fastq_ena(sra_id: str, outdir: str, outfiles: list,
                       method: str,
                       aspera_config: Optional[str],
                       cmdout: str, logi: TextIO):
    """
    Wrapper to download files using enaBrowserTools.
    There is a bug in the current version of the tool which breaks it in
    Python >= 3.11 so I have a local submodule version with two changes
    in the utils.py file to fix this (changed SafeConfigParser to
    ConfigParser).
    If sra_aspera is specified, use the Aspera connect tool to increase \
    (hopefully) download speed.

    Parameters
    ----------
    sra_id: str
        SRA run ID to download
    outdir: str
        Output directory for log files etc
    outfiles: list
        List of paths to single end, paired end 1, paired end 2 desired output
        files
    method: str
        Either ena to download directly from ENA or ena_aspera to use
        the Aspera connect tool.
    aspera_config: str
        Path to Aspera config file
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file
    """
    # Need to run the local version of enaDataGet - determine the relative
    # path to the directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = "/".join(script_dir.split("/")[:-2])
    ena_dir = f"{root_dir}/enaBrowserTools/python3"

    # Build the statement to get the files
    if method == 'ena':
        statement = f'{ena_dir}/enaDataGet -f fastq -d {outdir} {sra_id}'
    elif method == 'ena_aspera':
        # If Aspera is specified, check if the config file is in the right
        # place, else fail
        if not os.path.exists(aspera_config):
            SGErrors.no_aspera_config(aspera_config)
        statement = f'{ena_dir}/enaDataGet -f fastq -as {aspera_config} \
-d {outdir} {sra_id}'

    # Get the data
    SGSetup.run_log_command(statement, cmdout, outdir, 'ena', sra_id, logi)

    # Aspera makes a folder with additional log files - move it into outdir
    if method == 'ena_aspera':
        shutil.move(f'{outdir}/{sra_id}/logs',
                    f'{outdir}/fastqs/enaDataGet_logs')

    # Paths where ENA will have stored the files
    p1 = f'{outdir}/{sra_id}/{sra_id}_1.fastq.gz'
    p2 = f'{outdir}/{sra_id}/{sra_id}_2.fastq.gz'
    s = f'{outdir}/{sra_id}/{sra_id}.fastq.gz'

    outs, outp1, outp2 = outfiles

    # Move the files to the right place, delete empty directories
    if os.path.exists(p1) and os.path.exists(p2):
        statement = f'mv {p1} {outp1}'
        SGSetup.run_log_command(statement, cmdout)

        statement = f'mv {p2} {outp2}'
        SGSetup.run_log_command(statement, cmdout)

        statement = f'rm -rf {outdir}/{sra_id}'
        SGSetup.run_log_command(statement, cmdout)
        return ("paired")

    elif os.path.exists(s):
        statement = f'mv {s} {outs}'
        SGSetup.run_log_command(statement, cmdout)
        statement = f'rm -rf {outdir}/{sra_id}'
        SGSetup.run_log_command(statement, cmdout)
        return ("single")


def download_fastq(sra_id: str, method: str, outfiles: list, outdir: str,
                   aspera_config: Optional[str], cmdout: str, logi: TextIO):
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
    aspera_config: str
        Path to Aspera config file, only needed if using Aspera
    cmdout:
        Path to file to store bash commands
    logi: logging.Logger
        Path to open log file
    Returns
    -------
    None
    '''

    if method == 'sra':
        # Download using NCBI fasterq-dump
        typ = download_fastq_fasterqdump(sra_id, outdir, outfiles,
                                         cmdout, logi)
    elif method == 'ena' or method == 'ena_aspera':
        # Download using ENA enaDataGet
        typ = download_fastq_ena(sra_id, outdir, outfiles, method,
                                 aspera_config, cmdout, logi)
    else:
        SGErrors.bad_download_source(method)

    # Store whether the data is single or paired end in outdir/data_type.txt
    typ_out = open(f"{outdir}/data_type.txt", "w")
    typ_out.write(typ)
    typ_out.close()

    # Make empty placeholders for single / paired - whichever the data
    # is not
    if typ == 'paired':
        SGSetup.touch_only([outfiles[0]], cmdout, logi)
    else:
        SGSetup.touch_only([outfiles[1], outfiles[2]], cmdout, logi)


def trim_reads_trim_galore(in_s, in_p1, in_p2,
                           out_s, out_p1, out_p2,
                           outdir, prefix, typ, cmdout, logi):
    if typ == 'paired':
        statement = f"trim_galore --paired -o {outdir}/trimmed {in_p1} {in_p2}"
        SGSetup.run_log_command(statement, cmdout, outdir,
                                'trim_galore', prefix, logi)
        statement = f"mv {outdir}/trimmed/{prefix}.fq.1.gz_val_1.fq.gz \
                      {out_p1}"
        SGSetup.run_log_command(statement, cmdout)
        statement = f"mv {outdir}/trimmed/{prefix}.fq.2.gz_val_2.fq.gz \
                      {out_p2}"
        SGSetup.run_log_command(statement, cmdout)
        SGSetup.touch_only([out_s], cmdout, logi)
    else:
        statement = f"trim_galore -o {outdir}/trimmed {in_s}"
        SGSetup.run_log_command(statement, cmdout, outdir,
                                'trim_galore', prefix, logi)
        statement = f"mv {outdir}/trimmed/{prefix}_trimmed.fq.gz {out_s}"
        SGSetup.run_log_command(statement, cmdout)
        SGSetup.touch_only([out_p1, out_p2], cmdout, logi)


def trim_reads(infiles: list,
               outfiles: list,
               outdir: str,
               prefix: str,
               trimmer: str,
               typ: str,
               cmdout: str,
               logi: TextIO):
    in_s, in_p1, in_p2 = infiles
    out_s, out_p1, out_p2 = outfiles
    if trimmer == 'trim_galore':
        trim_reads_trim_galore(in_s, in_p1, in_p2,
                               out_s, out_p1, out_p2,
                               outdir, prefix, typ, cmdout, logi)
    else:
        SGErrors.trimmer_not_found(trimmer)
