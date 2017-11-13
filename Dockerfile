FROM andrewosh/binder-base

USER root

# Add Julia dependencies
RUN apt-get update && apt-get install samtools

USER main

