FROM continuumio/miniconda:4.7.12

RUN apt-get update && apt-get install -y gcc zlib1g-dev

WORKDIR /app
COPY environment.yml ./

RUN conda env create -f environment.yml
