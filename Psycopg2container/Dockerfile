FROM quay.io/biocontainers/biopython:1.79
LABEL maintainer="Tyler Peirce tyler.peirce@uwa.edu.au"
LABEL description="This is a docker file for running psycopg2 to update postgreSQL database"

# Set Environments
ENV PATH="/bin:${PATH:+:${PATH}}"

# Install app dependencies
RUN pip install psycopg2-binary
RUN pip install pandas
RUN pip install numpy

