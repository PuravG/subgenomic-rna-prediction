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

    Returns
    -------
    None
    '''
    # Links need an absolute path
    ap = os.path.abspath(infile)
    statement = f'pigz -c {ap} > {outfile}'
    SGSetup.run_log_command(cmdout, statement)


def zip_and_link(fastq_s: Optional[str],
                 fastq_p1: Optional[str],
                 fastq_p2: Optional[str],
                 outfiles: list,
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
                SGSetup.run_log_command(cmdout, statement)
            else:
                logi.info(f"FASTQ {fastq_s} is unzipped - zipping")
                # gzip the file
                zip_to_file(fastq_s, outfile_s, cmdout, logi)
            # Empty placeholders for paired end
            SGSetup.touch_only([outfile_p1, outfile_p2], cmdout, logi)
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
                zip_to_file(fastq_p1, outfile_p1, cmdout, logi)
                zip_to_file(fastq_p2, outfile_p2, cmdout, logi)
            # Empty placeholder for single end
            SGSetup.touch_only([outfile_s], cmdout, logi)
        else:
            # If the file path is not valid
            SGErrors.fastq_not_found([fastq_p1, fastq_p2])
    else:
        # If both file paths are None or False raise an error
        SGErrors.no_files_found()


def download_fastq_fasterqdump(sra_id: str, outdir: str, outfiles: list,
                               cmdout: str, logi: TextIO):
    # Download the files
    statement = f'fasterq-dump --split-files {sra_id} --outdir {outdir}'
    SGSetup.run_log_command(cmdout, statement, outdir, 'fasterq',
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
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = "/".join(script_dir.split("/")[:-2])
    ena_dir = f"{root_dir}/enaBrowserTools/python3"
    if method == 'ena':
        statement = f'{ena_dir}/enaDataGet -f fastq -d {outdir} {sra_id}'
    elif method == 'ena_aspera':
        statement = f'{ena_dir}/enaDataGet -f fastq -as {aspera_config} \
-d {outdir} {sra_id}'
    SGSetup.run_log_command(cmdout, statement, outdir, 'ena', sra_id, logi)
    if method == 'ena_aspera':
        shutil.move(f'{outdir}/{sra_id}/logs',
                    f'{outdir}/fastqs/enaDataGet_logs')
    p1 = f'{outdir}/{sra_id}/{sra_id}_1.fastq.gz'
    p2 = f'{outdir}/{sra_id}/{sra_id}_2.fastq.gz'
    s = f'{outdir}/{sra_id}/{sra_id}.fastq.gz'

    outs, outp1, outp2 = outfiles
    if os.path.exists(p1) and os.path.exists(p2):
        statement = f'mv {p1} {outp1}'
        SGSetup.run_log_command(cmdout, statement)

        statement = f'mv {p2} {outp2}'
        SGSetup.run_log_command(cmdout, statement)

        statement = f'rm -rf {outdir}/{sra_id}'
        SGSetup.run_log_command(cmdout, statement)
        return ("paired")

    elif os.path.exists(s):
        statement = f'mv {s} {outs}'
        SGSetup.run_log_command(cmdout, statement)
        statement = f'rm -rf {outdir}/{sra_id}'
        SGSetup.run_log_command(cmdout, statement)
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
    cmdout: str
        Path to command line log file
    logi: logger.
    Returns
    -------
    None
    '''
    typ_out = open(f"{outdir}/data_type.txt", "w")
    if method == 'sra':
        typ = download_fastq_fasterqdump(sra_id, outdir, outfiles,
                                         cmdout, logi)
    elif method == 'ena' or method == 'ena_aspera':
        typ = download_fastq_ena(sra_id, outdir, outfiles, method,
                                 aspera_config, cmdout, logi)
    else:
        SGErrors.bad_download_source(method)
    typ_out.write(typ)
    typ_out.close()
    if typ == 'paired':
        SGSetup.touch_only([outfiles[0]], cmdout, logi)
    else:
        SGSetup.touch_only([outfiles[1], outfiles[2]], cmdout, logi)
