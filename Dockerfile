FROM andrewosh/binder-base

USER root

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}

# Add samtools
RUN apt-get update && apt-get install samtools

USER main

RUN pip install -r requirements.txt

RUN gunzip data/chr2.fa.gz
