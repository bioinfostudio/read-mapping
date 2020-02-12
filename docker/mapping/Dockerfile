FROM asia.gcr.io/wise-cycling-267806/mapping-base
WORKDIR /home/bioinfo/mapping/
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.fastq
WORKDIR /home/bioinfo/mapping/bowtie_index/
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.1.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.2.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.3.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.4.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.fa
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.rev.1.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mm10.rev.2.bt2
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/mouse.mm10.genome
# setup the user directory
RUN chown -R bioinfo:bioinfo /home/bioinfo
RUN echo 'export PS1=" \W\$ "' >> /home/bioinfo/.bashrc
# become the bioinfo
USER bioinfo
ENV HOME /home/bioinfo
WORKDIR /home/bioinfo/
