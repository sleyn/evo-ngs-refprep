FROM ubuntu:20.04 AS build

LABEL author "Semion Leyn"
LABEL email "semen.leyn@gmail.com"
LABEL description "Preparation of the reference genome for the variant calling pipeline"

WORKDIR /root/

# Disable all questions from the apt
ENV DEBIAN_FRONTEND=noninteractive

# Create directory for binaries that should be copied to the new stage
RUN mkdir -p ./bin

# Install all dependences
RUN apt-get update \
    && apt-get -y install wget \
                          make \
                          gcc \
                          libbz2-dev \
                          libncurses5-dev \
                          # libncursesw5-dev \
                          csh \
                          libcurl4-openssl-dev \
                          liblzma-dev \
                          zlib1g-dev \
                          bwa \
                          # perl \
                          # cpanminus \
    && cp /usr/bin/bwa ./bin \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install bwa
# RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
#     && tar -xvjf bwa-0.7.17.tar.bz2 \
#     && cd bwa-0.7.17 \
#     && make \
#     && mkdir -p ../bin \
#     && mv bwa ../bin \
#     && cd .. \
#     && rm -rf bwa*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 \
    && tar -xvjf samtools-1.15.1.tar.bz2 \
    && cd samtools-1.15.1 \
    && make \
    && mv samtools ../bin \
    && cd ..

# Install picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.2/picard.jar \
    && mv picard.jar ./bin

# Install MUMmer
RUN wget https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz \
    && tar -xvzf MUMmer3.23.tar.gz \
    && cd MUMmer3.23 \
    && make \
    && cd .. \
    && mv MUMmer3.23 ./bin

# Configure main release
FROM ubuntu:20.04
WORKDIR /root/
COPY --from=build /root/bin/* /usr/local/bin/

# Add Mummer directory to PATH
ENV PATH=${PATH}:/usr/local/bin/MUMmer3.23

# Disable all questions from the apt
ENV DEBIAN_FRONTEND=noninteractive

# Install required files
RUN apt-get update \
    && apt-get -y install perl \
                          cpanminus \
                          python3 \
                          python3-pip \
                          libcurl4 \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Add JSON module to perl
RUN cpanm JSON 

# Add Pandas and Biopython
RUN pip3 install pandas biopython

ENTRYPOINT ["bash", "-l", "-c"]