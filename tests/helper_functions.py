import gzip
import hashlib
import subprocess


def matchRunTimeErrorSnakemake(snakemake_output, string):
    """
    Read the STDERR from a snakemake run and look for a specific error
    message.
    """
    for line in str(snakemake_output.stderr.decode('utf-8')).split("\n"):
        if string in line:
            return True
    return False


def compare_zipped(filepath1, filepath2):
    """
    Compare each line of two zipped files and check that all are identical
    """
    g1 = gzip.open(filepath1).readlines()
    g2 = gzip.open(filepath2).readlines()
    for l1, l2 in zip(g1, g2):
        assert l1 == l2
    return True


def hash_file(filepath):
    '''
    Retrieve the md5sum for a file
    '''
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
                -R --until %s" % (config, func)
    result = subprocess.run(statement, shell=True,
                            capture_output=True)
    return result
