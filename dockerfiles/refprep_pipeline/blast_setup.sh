#!/usr/bin/env bash

if [ ${TARGETPLATFORM} == "linux/amd64" ]; then \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && tar -xvzf ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && mv ncbi-blast-2.13.0+/bin/* ./bin
elif [ ${TARGETPLATFORM} == "linux/arm64" ]
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-arm-linux.tar.gz \
    && tar -xvzf ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && mv ncbi-blast-2.13.0+/bin/* ./bin
fi