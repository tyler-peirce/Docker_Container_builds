FROM pawsey/mpich-base:3.1.4_ubuntu20.04
LABEL maintainer="Tyler Peirce tyler.peirce@uwa.edu.au"
LABEL description="This is a docker file for running mitogenome annotation program EMMA and coding language Julia"

ENV HOME=/opt

# set environments
ENV ROCM_PATH=/opt/rocm
ENV MPICH_HOME="/opt/mpich"
ENV PATH="/opt/julia-1.10.5/bin:bin:/opt/Emma:${PATH:+:${PATH}}"
ENV MANPATH="/usr/share/man/man1:${MANPATH:+:${MANPATH}}"
ENV JULIA_DEPOT_PATH="/opt/julia-1.10.5/bin" 

# install app dependencies
RUN apt-get update && apt-get install -y \
        hmmer \
        infernal infernal-doc \
        git \
        lua5.3 \
        curl

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz
RUN tar zxvf julia-1.10.5-linux-x86_64.tar.gz
RUN mv julia-1.10.5 /opt/.

RUN git clone https://github.com/ian-small/Emma /opt/Emma

# Set up EMMA as packages
RUN /opt/julia-1.10.5/bin/julia -e 'using Pkg; \
        Pkg.develop(PackageSpec(path="/opt/Emma")); \
        Pkg.instantiate()'
        
RUN /opt/julia-1.10.5/bin/julia -e 'using Pkg; \
        Pkg.add("FASTX"); \
        Pkg.add("BioSequences")'

# Precompile Julia packages (optional, can speed up runtime)
RUN julia -e 'using Pkg; Pkg.precompile()'


COPY extract_proteins.jl /opt/.

RUN chmod +x /opt/*.jl



