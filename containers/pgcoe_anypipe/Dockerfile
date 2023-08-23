FROM continuumio/miniconda3:22.11.1

# switch to mamba & install conda packages
WORKDIR /root/
COPY environment.yml .

ENV PATH /opt/conda/envs/env/bin:$PATH
RUN conda config --add channels bioconda
RUN echo $(conda config --show channels)
RUN conda env update -n root -f environment.yml && conda clean -a

#BCFtools / HTSlib / Samtools
RUN conda install -c bioconda -c conda-forge -y htslib'>=1.13,<1.14.0a0'
RUN conda install -c bioconda -c conda-forge -y bcftools=1.13
RUN conda install -c bioconda -c conda-forge -y samtools=1.11
RUN conda install -c bioconda -c conda-forge -y bedtools


ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"
