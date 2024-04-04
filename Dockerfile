FROM ubuntu:20.04
# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Install dependencies, including liblzma and other previously mentioned packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    less \
    unzip \
    gzip \
    make \
    zstd \
    bash \
    git \
    python3 \
    python3-pip \
    python3-dev \
    samtools \
    bcftools \
    gcc \
    g++ \
    libc-dev \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz --no-check-certificate \
    && tar -xvzf BBMap_38.90.tar.gz \
    && rm BBMap_38.90.tar.gz \
    && mv bbmap /opt/
# Update PATH for BBTools
ENV PATH="/opt/bbmap:${PATH}"
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_23.10.0-1-Linux-aarch64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

RUN conda config --add channels conda-forge && \
    conda update -n base --all && \
    conda install -y -n base \
    mamba \
    libmamba \
    libmambapy \
    conda-libmamba-solver

# COPY requirements.yml /tmp/requirements.yml
RUN mamba install -y -n base -c bioconda -c conda-forge pandas biopython pandas numpy nextflow openpyxl scipy scikit-learn

COPY ./main.nf /scratch/main.nf
COPY ./nextflow.config /scratch/nextflow.config

# make sure shells are bash
CMD ["/bin/bash"]