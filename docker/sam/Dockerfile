FROM asia.gcr.io/wise-cycling-267806/mapping

ADD https://storage.googleapis.com/bioinfostudio/mapping/data/Oct4.sam /home/bioinfo/mapping/

USER root
RUN chown -R bioinfo:bioinfo /home/bioinfo
USER bioinfo
WORKDIR /home/bioinfo/mapping/
