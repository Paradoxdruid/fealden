# Two-stage Dockerfile with binary builder and compiled final stage

FROM ubuntu:22.04 AS builder

ENV HOME /root

SHELL ["/bin/bash", "-c"]

# Install build system necessities
RUN apt-get update && apt-get -y --no-install-recommends install \
    curl=7.81.0-1ubuntu1.10 \
    python3.10=3.10.6-1~22.04.2ubuntu1 \
    python3-pip=22.0.2+dfsg-1ubuntu0.2 \
    libpython3.10-dev=3.10.6-1~22.04.2ubuntu1 \
    wget=1.21.2-2ubuntu1 \
    build-essential=12.9ubuntu3 \
    swig=4.0.2-1ubuntu1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies to the local user directory (eg. /root/.local)
RUN pip install --user python-dotenv==1.0.0 pytest==7.2.2

WORKDIR ${HOME}

# Install UNAfold
RUN wget --no-check-certificate --quiet \
    https://rnaspace.sourceforge.net/software/unafold-3.8.tar.gz && \
    tar xzf ./unafold-3.8.tar.gz && \
    cd ./unafold-3.8 && \
    ./configure --prefix=/root/unafold && \
    make && \
    make install

# Install mfold
RUN wget --no-check-certificate --quiet \
    http://www.unafold.org/download/mfold-3.6.tar.gz && \
    tar xzf ./mfold-3.6.tar.gz && \
    cd ./mfold-3.6 && \
    ./configure --prefix=/root/mfold && \
    make && \
    make install

# Working stage
FROM ubuntu:22.04

ENV HOME /root

SHELL ["/bin/bash", "-c"]

# install build system necessities
RUN apt-get update && apt-get -y --no-install-recommends install \
    python3.10=3.10.6-1~22.04.2ubuntu1 \
    python3-pip=22.0.2+dfsg-1ubuntu0.2 \
    wget=1.21.2-2ubuntu1 \
    curl=7.81.0-1ubuntu1.10 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Binary dependencies
COPY --from=builder /root/mfold /root/mfold
COPY --from=builder /root/unafold /root/unafold

# Local pip installs
COPY --from=builder /root/.local /root/.local
ENV PATH=/root/.local:$PATH

RUN ln -s /usr/bin/python3 /usr/bin/python

# Fealden environment variables
ENV FEALDEN_BACKEND=mfold
ENV HYBRID_SS_MIN=/root/unafold/bin/hybrid-ss-min
ENV SIR_GRAPH=/root/mfold/bin/sir_graph

# Get latest fealden
ARG FORCE_UPDATE=no
RUN cd ${HOME} && \
    wget --no-check-certificate \
    https://github.com/Paradoxdruid/fealden/archive/refs/heads/master.tar.gz \
    -O "fealden-master.tar.gz" && \
    mkdir ./fealden && \
    tar -xzvf ./"fealden-master.tar.gz" -C ./fealden --strip-components 1 && \
    rm ./"fealden-master.tar.gz"

WORKDIR /root/fealden
