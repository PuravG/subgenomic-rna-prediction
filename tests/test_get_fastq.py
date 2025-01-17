import pytest
import os
import shutil
import helper_functions


@pytest.mark.parametrize(
    "test_file, config_file, stem",
    [("se_zip.fq.gz",
      "config_se_zipped.yaml",
      "fastqs_se_zip"),
     ("se_unzip.fq.gz",
      "config_se_unzipped.yaml",
      "fastqs_se_unzip")])
def test_getfastq_se(test_file, config_file, stem):
    """Test Snakefile for user provided single end zipped and unzipped \
       FASTQ files"""
    # Set up file paths
    test_output_dir = f"tests/test_output_{stem}"
    test_output = f"{test_output_dir}/fastqs/{test_file}"
    expected_dir = "tests/expected_output/fastq_tests"
    expected_output = f"{expected_dir}/SRR7660720_small.fastq.gz"
    config_path = "tests/config_files/fastq_tests/%s" % config_file

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0  # Check if snakemake ran successfully

    # Check the output files are identical to the expected outputs
    assert helper_functions.compare_zipped(expected_output, test_output)

    # Check the paired end placeholders exist
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.1.gz"))
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.2.gz"))
    shutil.rmtree(test_output_dir)


@pytest.mark.parametrize(
    "test_files, config_file, stem",
    [(["pe_zip.fq.1.gz",
       "pe_zip.fq.2.gz"],
      "config_pe_zipped.yaml",
      "fastqs_pe_zip"),
     (["pe_unzip.fq.1.gz",
       "pe_unzip.fq.2.gz"],
      "config_pe_unzipped.yaml",
      "fastqs_pe_unzip")])
def test_getfastq_pe(test_files, config_file, stem):
    """Test Snakefile for user provided paired end zipped and unzipped \
       FASTQ files"""
    # Set up file paths
    test_output_dir = f"tests/test_output_{stem}"
    test_file_1, test_file_2 = test_files
    test_output_1 = f"{test_output_dir}/fastqs/{test_file_1}"
    test_output_2 = f"{test_output_dir}/fastqs/{test_file_2}"
    expected_dir = "tests/expected_output/fastq_tests"
    expected_output_1 = f"{expected_dir}/SRR18609750_1_small.fastq.gz"
    expected_output_2 = f"{expected_dir}/SRR18609750_2_small.fastq.gz"
    config_path = "tests/config_files/fastq_tests/%s" % config_file

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0

    # Check the output files are identical to the expected outputs
    assert helper_functions.compare_zipped(expected_output_1, test_output_1)
    assert helper_functions.compare_zipped(expected_output_2, test_output_2)

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
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "At least one of your FASTQ files does not exist")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


def test_getfastq_allmissing():
    """Test Snakefile with missing input fastq files - should raise
       a runtime error -
       FASTQ input specified but no fastq files provided"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/config_missing_all.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "FASTQ input specified but no fastq files provided")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


def test_getfastq_bad_start_point():
    """Test Snakefile with wrong start point - not raw_fastq,
    trimmed fastq, sra_id or bam. Should raise a Runtime error"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/config_bad_start_point.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "The start_point parameter should be")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


@pytest.mark.parametrize(
    "config_file, stem, suffix",
    [('config_pe_sra', 'fastqs_pe_sra', 'sra'),
     ('config_pe_ena', 'fastqs_pe_ena', 'ena'),
     ('config_pe_ena_aspera', 'fastqs_pe_ena_aspera', 'ena')])
def test_getfastq_online_paired(config_file, stem, suffix):
    # Set up file paths
    expected_dir = "tests/expected_output/fastq_tests"
    expected_output_1 = f"{expected_dir}/ERR5715266_{suffix}.fq.1.gz"
    expected_output_2 = f"{expected_dir}/ERR5715266_{suffix}.fq.2.gz"
    config_path = "tests/config_files/fastq_tests/%s.yaml" % config_file
    test_output_dir = f"tests/test_output_{stem}"
    test_output_1 = f"{test_output_dir}/fastqs/ERR5715266.fq.1.gz"
    test_output_2 = f"{test_output_dir}/fastqs/ERR5715266.fq.2.gz"
    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")
    assert result.returncode == 0

    # Check the output is correct
    assert helper_functions.compare_zipped(expected_output_1, test_output_1)
    assert helper_functions.compare_zipped(expected_output_2, test_output_2)
    shutil.rmtree(test_output_dir)


@pytest.mark.parametrize(
    "config_file, stem, suffix",
    [('config_se_sra', 'fastqs_se_sra', 'sra'),
     ('config_se_ena', 'fastqs_se_ena', 'ena'),
     ('config_se_ena_aspera', 'fastqs_se_ena_aspera', 'ena')])
def test_getfastq_online_single(config_file, stem, suffix):
    # Set up file paths
    expected_dir = "tests/expected_output/fastq_tests"
    expected_output = f"{expected_dir}/SRR28628208_{suffix}.fq.gz"
    config_path = "tests/config_files/fastq_tests/%s.yaml" % config_file
    test_output_dir = f"tests/test_output_{stem}"
    test_output = f"{test_output_dir}/fastqs/SRR28628208.fq.gz"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check the output is correct
    assert result.returncode == 0

    assert helper_functions.compare_zipped(expected_output, test_output)
    shutil.rmtree(test_output_dir)


def test_getfastq_bad_source():
    """Test Snakefile with wrong source - not sra, ena, ena_aspera"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/config_bad_source.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "Download method")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


@pytest.mark.parametrize("config_file",
                         ["config_bad_sra",
                          "config_bad_ena"])
def test_getfastq_bad_ID(config_file):
    """Test Snakefile with invalid SRA ID and sra"""
    # Set up file paths
    config_path = "tests/config_files/fastq_tests/%s.yaml" % config_file

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "run failed with ID ")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


def test_getfastq_no_aspera_config():
    """Test Snakefile with aspera specified but no aspera config file"""
    # Set up file paths
    config_dir = 'tests/config_files/fastq_tests'
    config_path = f"{config_dir}/config_se_ena_aspera_noconfig.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_fastq")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "Aspera config file not found at")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


def test_getfastq_skip():
    """Test Snakefile with bam as the start_point - should just make empty
    placeholder files"""
    config_dir = 'tests/config_files/fastq_tests'
    config_path = f"{config_dir}/config_bam_start.yaml"
    test_output_dir = "tests/test_output_fastqs_bam_start"
    test_output_p1 = f"{test_output_dir}/fastqs/bam_start.fq.1.gz"
    test_output_p2 = f"{test_output_dir}/fastqs/bam_start.fq.2.gz"
    test_output_s = f"{test_output_dir}/fastqs/bam_start.fq.gz"

    result = helper_functions.run_snakemake(config_path, "get_fastq")
    assert result.returncode == 0, result.returncode
    assert os.path.getsize(test_output_p1) == 0
    assert os.path.getsize(test_output_p2) == 0
    assert os.path.getsize(test_output_s) == 0

    result = helper_functions.run_snakemake(config_path, "trim_fastq")
    assert result.returncode == 0, result.returncode
    assert os.path.getsize(test_output_p1) == 0
    assert os.path.getsize(test_output_p2) == 0
    assert os.path.getsize(test_output_s) == 0
    shutil.rmtree('tests/test_output_fastqs_bam_start')


def testTrimPE():
    """Test trimming for paired end files"""
    config_dir = 'tests/config_files/fastq_tests'
    config_path = f"{config_dir}/config_pe_trim.yaml"
    test_output_dir = "tests/test_output_fastqs_pe_trim"

    expected_dir = "tests/expected_output/fastq_tests"
    expected_output_1 = f"{expected_dir}/SRR18609750_1_trimmed.fastq.gz"
    expected_output_2 = f"{expected_dir}/SRR18609750_2_trimmed.fastq.gz"

    test_output_1 = f"{test_output_dir}/trimmed/pe_trim.fq.1.gz"
    test_output_2 = f"{test_output_dir}/trimmed/pe_trim.fq.2.gz"

    result = helper_functions.run_snakemake(config_path, "trim_fastq")
    assert result.returncode == 0, result.returncode

    assert helper_functions.hash_file(
        expected_output_1) == helper_functions.hash_file(test_output_1)
    assert helper_functions.hash_file(
        expected_output_2) == helper_functions.hash_file(test_output_2)
    shutil.rmtree('tests/test_output_fastqs_pe_trim')


def testTrimSE():
    """Test trimming for single end files"""
    config_dir = 'tests/config_files/fastq_tests'
    config_path = f"{config_dir}/config_se_trim.yaml"
    test_output_dir = "tests/test_output_fastqs_se_trim"

    expected_dir = "tests/expected_output/fastq_tests"
    expected_output = f"{expected_dir}/SRR7660720_trimmed.fastq.gz"

    test_output = f"{test_output_dir}/trimmed/se_trim.fq.gz"

    result = helper_functions.run_snakemake(config_path, "trim_fastq")
    assert result.returncode == 0, result.returncode

    assert helper_functions.hash_file(
        expected_output) == helper_functions.hash_file(test_output)
    shutil.rmtree('tests/test_output_fastqs_se_trim')


def testTrimFail():
    """Test trimming with a trimmer that is not implemented - should fail"""
    config_dir = 'tests/config_files/fastq_tests'
    config_path = f"{config_dir}/config_bad_trim.yaml"

    result = helper_functions.run_snakemake(config_path, "trim_fastq")
    assert result.returncode != 0
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "Trimming tool ")
    shutil.rmtree('tests/test_output_fastqs_bad_trim')


def test_getfastq_pe_trimmed_existing():
    """Test Snakefile for user provided trimmed paired end zipped and \
       unzipped FASTQ files"""
    # Set up file paths
    test_output_dir = "tests/test_output_fastqs_pe_trimmed_existing"
    test_output_1 = f"{test_output_dir}/trimmed/pe_trimmed_existing.fq.1.gz"
    test_output_2 = f"{test_output_dir}/trimmed/pe_trimmed_existing.fq.2.gz"

    expected_dir = "tests/expected_output/fastq_tests"
    expected_output_1 = f"{expected_dir}/SRR18609750_1_small.fastq.gz"
    expected_output_2 = f"{expected_dir}/SRR18609750_2_small.fastq.gz"

    config_dir = "tests/config_files/fastq_tests"
    config_path = f"{config_dir}/config_pe_trimmed_existing.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "trim_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0

    # Check the output files are identical to the expected outputs
    assert helper_functions.compare_zipped(expected_output_1, test_output_1)
    assert helper_functions.compare_zipped(expected_output_2, test_output_2)

    # Check the single end placeholder exists
    assert os.path.exists(test_output_1.replace(".fq.1.gz", ".fq.gz"))

    # Remove test data
    shutil.rmtree(test_output_dir)


def test_getfastq_se_trimmed_existing():
    """Test Snakefile for user provided trimmed paired end zipped and \
       unzipped FASTQ files"""
    # Set up file paths
    test_output_dir = "tests/test_output_fastqs_se_trimmed_existing"
    test_output = f"{test_output_dir}/trimmed/se_trimmed_existing.fq.gz"

    expected_dir = "tests/expected_output/fastq_tests"
    expected_output = f"{expected_dir}/SRR7660720_small.fastq.gz"

    config_dir = "tests/config_files/fastq_tests"
    config_path = f"{config_dir}/config_se_trimmed_existing.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "trim_fastq")

    # Check if snakemake ran successfully
    assert result.returncode == 0

    # Check the output files are identical to the expected outputs
    assert helper_functions.compare_zipped(expected_output, test_output)

    # Check the single end placeholder exists
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.1.gz"))
    assert os.path.exists(test_output.replace(".fq.gz", ".fq.2.gz"))

    # Remove test data
    shutil.rmtree(test_output_dir)
