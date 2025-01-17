#!/usr/bin/env python3
import os
import SGErrors
import SGSetup
from typing import Optional
from typing import TextIO
import shutil


def link_genome(genome_path: str,
                genome_name: str,
                outdir: str,
                cmdout: str,
                logi: TextIO):
    if not genome_path:
        SGErrors.no_genome_spec()
    if not os.path.exists(genome_path):
        SGErrors.no_genome_file(genome_path)
