FROM andrewosh/binder-base

USER root

# Add samtools
RUN apt-get update && apt-get install samtools

USER main

RUN pip install -r requirements.txt



