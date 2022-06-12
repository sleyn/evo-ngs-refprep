#!/usr/bin/env bash

if [ ${TARGETPLATFORM} == "linux/amd64" ]; then
    TARFILE="ncbi-blast-2.13.0+-x64-linux.tar.gz"
elif [ ${TARGETPLATFORM} == "linux/arm64" ]; then 
    TARFILE="ncbi-blast-2.13.0+-x64-arm-linux.tar.gz"
fi

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/${TARFILE}
tar -xvzf ${TARFILE}
mv ncbi-blast-2.13.0+/bin/* /usr/local/bin