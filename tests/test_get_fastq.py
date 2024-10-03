import pytest
import subprocess
import os
import shutil
import hashlib
import gzip


def matchRunTimeErrorSnakemake(snakemake_output, string):
    for line in str(snakemake_output.stderr.decode('utf-8')).split("\n"):
        if string in line:
            return True
    return False


def compare_zipped(filepath1, filepath2):
    g1 = gzip.open(filepath1).readlines()
    g2 = gzip.open(filepath2).readlines()
    for l1, l2 in zip(g1, g2):
        assert l1 == l2
    return True


def hash_file(filepath):
    hasher = hashlib.md5()  # You can use other algorithms like md5 or sha1
    with open(filepath, 'rb') as f:
        content = f.read()
        hasher.update(content)
    # Return the hexadecimal digest of the hash
    return hasher.hexdigest()


def run_snakemake(config, func):
    """Helper function to run snakemake."""
    statement = "snakemake \
                --snakefile subgenomic-rna-prediction/pipeline/Snakefile \
                --configfile %s \
                -R %s" % (config, func)
    result = subprocess.run(statement, shell=True,
                            capture_output=True)
    return result


@pytest.mark.parametrize(
    "test_file, exp_file, config_file, stem",
    [("se_zip.fq.gz",
      "SRR7660720_small.fastq.gz",
      "config_se_zipped.yaml",
      "fastqs_se_zip"),
     ("se_unzip.fq.gz",
      "SRR7660720_small.fastq.gz",
      "config_se_unzipped.yaml",
      "fastqs_se_unzip")])
def test_getfastq_se(test_file, exp_file, config_file, stem):
    """Test Snakefile for user provided single end zipped and unzipped \
       FASTQ files"""
    # Set up file paths
    test_output_dir = "tests/test_output_%s" % stem
    test_output = "%s/fastqs/%s" % (test_output_dir, test_file)
    expected_output = "tests/expected_output/fastq_tests/%s" % exp_file
    config_path = "tests/config_files/fastq_tests/%s" % config_file

    # Run the pipeline
    result = run_snakemake(config_path, "get_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0  # Check if snakemake ran successfully

    # Check the output files are identical to the expected outputs
    assert compare_zipped(expected_output, test_output)

    # Check the paired end placeholders exist
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.1.gz"))
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.2.gz"))
    shutil.rmtree(test_output_dir)


@pytest.mark.parametrize(
    "test_files, exp_files, config_file, stem",
    [(["pe_zip.fq.1.gz",
       "pe_zip.fq.2.gz"],
      ["SRR18609750_1_small.fastq.gz",
       "SRR18609750_2_small.fastq.gz"],
      "config_pe_zipped.yaml",
      "fastqs_pe_zip"),
     (["pe_unzip.fq.1.gz",
       "pe_unzip.fq.2.gz"],
      ["SRR18609750_1_small.fastq.gz",
       "SRR18609750_2_small.fastq.gz"],
      "config_pe_unzipped.yaml",
      "fastqs_pe_unzip")])
def test_getfastq_pe(test_files, exp_files, config_file, stem):
    """Test Snakefile for user provided paired end zipped and unzipped \
       FASTQ files"""
    # Set up file paths
    test_output_dir = "tests/test_output_%s" % stem
    test_output_1 = "%s/fastqs/%s" % (test_output_dir, test_files[0])
    test_output_2 = "%s/fastqs/%s" % (test_output_dir, test_files[1])
    expected_output_1 = "tests/expected_output/fastq_tests/%s" % exp_files[0]
    expected_output_2 = "tests/expected_output/fastq_tests/%s" % exp_files[1]
    config_path = "tests/config_files/fastq_tests/%s" % config_file

    # Run the pipeline
    result = run_snakemake(config_path, "get_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0

    # Check the output files are identical to the expected outputs
    assert compare_zipped(expected_output_1, test_output_1)
    assert compare_zipped(expected_output_2, test_output_2)

    # Check the single end placeholder exists
    assert os.path.exists(test_output_1.replace(".fq.1.gz", ".fq.gz"))

    # Remove test data
    shutil.rmtree(test_output_dir)


@pytest.mark.parametrize("config_file",
                         ["config_se_missing.yaml",
                          "config_pe_missing1.yaml",
                          "config_pe_missing2.yaml"])
def test_getfastq_missing(config_file):
    """Test Snakefile with missing input fastq files - should raise \
       a runtime error - One of your FASTQ files does not exist"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/%s" % config_file

    # Run the pipeline
    result = run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert matchRunTimeErrorSnakemake(
        result, "At least one of your FASTQ files does not exist")
    assert result.returncode == 1, result.returncode


def test_getfastq_allmissing():
    """Test Snakefile with missing input fastq files - should raise
       a runtime error -
       FASTQ input specified but no fastq files provided"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/config_missing_all.yaml"

    # Run the pipeline
    result = run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert matchRunTimeErrorSnakemake(
        result, "FASTQ input specified but no fastq files provided")
    assert result.returncode == 1, result.returncode


def test_getfastq_bad_start_point():
    """Test Snakefile with wrong start point - not raw_fastq,
    trimmed fastq, sra_id or bam"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/config_bad_start_point.yaml"

    # Run the pipeline
    result = run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert matchRunTimeErrorSnakemake(
        result, "The start_point parameter should be")
    assert result.returncode == 1, result.returncode
