FROM continuumio/miniconda:4.7.12

RUN apt-get update && apt-get install -y gcc zlib1g-dev

WORKDIR /app
COPY environment.yml ./

RUN conda env create -f environment.yml

RUN conda config --add channels bioconda && conda install trimmomatic=0.32 star=2.5

# copy in ribodiff dep
COPY ./RiboDiff /app/RiboDiff
RUN cd /app/RiboDiff && pip install -r requirements.txt && python setup.py build test install && cd /app

# copy in the rest of the source
COPY . ./

# RUN mkdir /usr/share/man/man1/ && apt-get install -y openjdk-11-jre-headless && \
#     curl -o /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip && \
#     unzip Trimmomatic-0.32.zip && \
#     cp /tmp/Trimmomatic-0.32/trimmomatic-0.32.jar /app

# CMD ["bash", "-c", "/app/entrypoint.sh"]s
COPY .bashrc /root/.bashrc
CMD /app/entrypoint.sh
