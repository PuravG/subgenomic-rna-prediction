import pytest
import os
import shutil
import helper_functions

config_dir = "tests/config_files/genome_tests"


def test_getgenome_no_fasta():
    """Should raise an error if virus_start_point is local but no
    genome FASTA file is specific"""
    # Set up file paths

    config_path = f"{config_dir}/config_getgenome_no_genome_path.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_genome")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "Genome FASTA file must be specified")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')


def test_getgenome_bad_fasta():
    """Should raise an error if virus_start_point is local but no
    genome FASTA file is specific"""
    # Set up file paths

    config_path = f"{config_dir}/config_getgenome_bad_genome_path.yaml"

    # Run the pipeline
    result = helper_functions.run_snakemake(config_path, "get_genome")

    # Check error is correct
    assert helper_functions.matchRunTimeErrorSnakemake(
        result, "Genome FASTA file not found at")
    assert result.returncode == 1, result.returncode
    shutil.rmtree('placeholder')