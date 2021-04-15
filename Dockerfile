FROM ubuntu:latest

MAINTAINER Todd N. Wylie

LABEL description = "Stack for BVI RNA-seq quantification."

ARG DEBIAN_FRONTEND=noninteractive

# STACK:
#  + FASTQC
#  + Kallisto
#  + MultiQC
#  + Pandas (Python)
#  + Python3
#  + Snakemake
#  + graphviz

# Update core applications and utilities.
RUN \
apt-get update && \
apt-get install -y \
apt-utils \
build-essential  \
make \
curl \
git \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
autoconf \
cmake \
unzip \
default-jre \
tar \
wget \
zsh \
python3.7-dev \
python3-pip \
python3-pandas \
graphviz \
zile

# Kallisto ####################################################################

RUN mkdir /usr/local/bin/installKallisto  && \
cd /usr/local/bin/installKallisto && \
wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
tar xvfz kallisto_linux-v0.46.1.tar.gz && \
cd kallisto/ && \
cp kallisto /usr/bin/

# FASTQC ######################################################################

RUN mkdir /usr/local/bin/installFASTQC && \\
cd /usr/local/bin/installFASTQC && \
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
unzip fastqc_v0.11.9.zip && \
cd FastQC && \
chmod +x fastqc && \
ln -s /usr/local/bin/installFASTQC/FastQC/fastqc /usr/bin/fastqc

# MultiQC #####################################################################

RUN python3.7 -m pip install --upgrade pip && \\
pip3 install --ignore-installed multiqc

# Snakemake ###################################################################

RUN python3.7 -m pip install snakemake==5.25.0

# Custom Code #################################################################

COPY bin/* /usr/bin/

# __END__
