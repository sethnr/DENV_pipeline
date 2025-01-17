from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from pathag_pipeline import __version__, _program

setup(name='pathag_pipeline',
      version=__version__,
      packages=find_packages() + ["pathag_pipeline.refs", "pathag_pipeline.snakemake"],
      scripts=[],
      description='Bioinformatic pipeline to generate reads and consensus sequences for any given viral genome',
      package_data={"pathag_pipeline.refs":["**"],
                    "pathag_pipeline.snakemake":["**"]},
      package_dir={
          'pathag_pipeline': 'pathag_pipeline',
          'pathag_pipeline.refs': './refs/',
          'pathag_pipeline.snakemake':"./snakemake",
      },
      install_requires=["biopython>=1.70"],
      url='https://github.com/ViralVerity/DENV_pipeline',
      author='Verity Hill',
      author_email='verity.hill@yale.edu',
      entry_points="""
          [console_scripts]
          {program} = pathag_pipeline.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
