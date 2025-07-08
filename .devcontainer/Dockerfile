# Base image with Ubuntu
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Fix locale warning
RUN apt-get update && apt-get install -y locales && \
    locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Install dependencies and tools
RUN apt-get update && \
    apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libtbb12 \
    python3 \
    python3-pip \
    default-jre \
    fastqc \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install fastp (precompiled binary)
#RUN wget http://opengene.org/fastp/fastp -O /usr/local/bin/fastp && \
    #chmod a+x /usr/local/bin/fastp


# Install STAR (precompiled binary)
#RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz && \
    #tar -xvzf 2.7.10a.tar.gz && \
    #cd STAR-2.7.10a/source && \
    #make STAR && \
    #cp STAR /usr/local/bin/ && \
    #cd ../.. && rm -rf STAR-2.7.10a 2.7.10a.tar.gz
   

# Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar -xvzf salmon-1.10.0_linux_x86_64.tar.gz && \
    mv salmon-latest_linux_x86_64/bin/salmon /usr/local/bin && \
    rm -rf salmon-latest_linux_x86_64*

# Set default shell
SHELL ["/bin/bash", "-c"]



