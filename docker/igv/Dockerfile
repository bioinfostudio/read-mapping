FROM asia.gcr.io/wise-cycling-267806/mapping

USER root
WORKDIR /home/bioinfo/mapping/

RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.sam
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.bam
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.sorted.bam
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.sorted.bam.bai
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.bedgraph
RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.bw

RUN chown -R bioinfo:bioinfo /home/bioinfo
USER bioinfo
