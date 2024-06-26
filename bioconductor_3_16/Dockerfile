FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e 'BiocManager::install("dada2", version="3.16", update=TRUE, ask=FALSE)'
RUN R -e 'BiocManager::install("decontam", version="3.16", update=TRUE, ask=FALSE)'

RUN R -e 'BiocManager::install("phyloseq", version="3.16", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("devtools", update=TRUE, ask=FALSE)'

RUN R -e 'BiocManager::install("ShortRead", version="3.16", update=TRUE, ask=FALSE)'
RUN R -e 'BiocManager::install("ANCOMBC", version="3.16", update=TRUE, ask=FALSE)'

RUN R -e 'devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")'
RUN R -e 'devtools::install_github("nuriamw/micro4all")'

RUN R -e 'BiocManager::install("Biostrings", version="3.16", update=TRUE, ask=FALSE)'
RUN R -e 'BiocManager::install("annotate", version="3.16", update=TRUE, ask=FALSE)'

RUN R -e 'BiocManager::install("edgeR", version="3.16", update=TRUE, ask=FALSE)'
RUN R -e 'BiocManager::install("microbiome", version="3.16", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("GUniFrac", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("phangorn", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("vegan", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("pheatmap", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("colorspace", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("ggplot2", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("tidyverse", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("seqinr", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("gdata", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("randomcoloR", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("scales", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("agricolae", update=TRUE, ask=FALSE)'

RUN R -e 'install.packages("svglite", update=TRUE, ask=FALSE)'

RUN apt-get update && \
    apt-get -y install cutadapt ncbi-blast+ mafft fasttree

# Exponer el puerto 8787
EXPOSE 8787

WORKDIR /opt
RUN git clone https://github.com/knono/figaro_py3.10.git

WORKDIR /opt/figaro_py3.10
RUN pip3 install --upgrade pip && \
    pip3 install -r requirements.txt

ENV PYTHONUNBUFFERED=1
