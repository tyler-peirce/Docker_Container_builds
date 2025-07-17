FROM pawsey/hpc-python:2021.09-hdf5mpi
LABEL maintainer="Tyler Peirce tyler.peirce@uwa.edu.au"
LABEL description="This is a docker file for running LCA using taxonkit, the environment TAXONKIT_DB can be definied to change the ~/.taxonkit location"

# Set Environments
ENV PATH="/bin:${PATH:+:${PATH}}"

# Install app dependencies
RUN apt-get update && apt-get install -y \
        procps

# Install TaxonKit
RUN wget https://github.com/shenwei356/taxonkit/releases/download/v0.17.0/taxonkit_linux_amd64.tar.gz
RUN tar -zxvf taxonkit_linux_amd64.tar.gz

# Copy the LCA script into the container
COPY computeLCA.py /.