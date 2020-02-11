FROM asia.gcr.io/wise-cycling-267806/mapping-base

ADD https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.fastq /home/bioinfo/mapping/

ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.1.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.2.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.3.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.4.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.fa /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.rev.1.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.rev.2.bt2 /home/bioinfo/mapping/bowtie_index/
ADD https://storage.googleapis.com/bioinfostudio/mapping/data/mouse.mm10.genome /home/bioinfo/mapping/bowtie_index/

# setup the user directory
RUN chown -R bioinfo:bioinfo /home/bioinfo
RUN echo 'export PS1=" \W\$ "' >> /home/bioinfo/.bashrc
# become the bioinfo
USER bioinfo
ENV HOME /home/bioinfo
WORKDIR /home/bioinfo/
