FROM asia.gcr.io/wise-cycling-267806/mapping

USER root
WORKDIR /home/bioinfo/mapping/

RUN wget https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.sam

RUN chown -R bioinfo:bioinfo /home/bioinfo
USER bioinfo
