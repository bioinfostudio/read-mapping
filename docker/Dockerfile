FROM biocontainers/biocontainers:vdebian-buster-backports_cv1

USER root
RUN apt-get update --fix-missing
RUN apt-get install -t buster-backports -y \
    curl \
    git \
    build-essential \
    nano \
    man \
    less \
    bowtie2 \
    samtools \
    bedtools

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph /usr/local/bin 
RUN chmod +x /usr/local/bin/bigWigToBedGraph
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# create a non-root user
RUN useradd -ms /bin/bash bioinfo 
