FROM biocontainers/samtools:1.3.1

#RUN apt-get update
#RUN pip install --no-cache-dir notebook==5.*


ENV NB_USER tt
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
    
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN gunzip data/chr2.fa.gz
