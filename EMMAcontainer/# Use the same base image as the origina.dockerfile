# Use the same base image as the original (you may need to adjust this)
FROM ubuntu:20.04

# install app dependencies
RUN apt-get update && apt-get install -y \
        hmmer \
        infernal infernal-doc \
        git \
        lua5.3 \
        curl \
        wget \
        && rm -rf /var/lib/apt/lists/*

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz \
    && tar -xzf julia-1.10.5-linux-x86_64.tar.gz \
    && mv julia-1.10.5 /opt/julia-1.10.5 \
    && rm julia-1.10.5-linux-x86_64.tar.gz

# Add Julia to PATH
ENV PATH="/opt/julia-1.10.5/bin:$PATH"
# Set Julia load path to include Emma
ENV JULIA_LOAD_PATH="/opt/Emma:/opt/Emma/src:$JULIA_LOAD_PATH"
# Make sure the Julia environment is activated by default
ENV JULIA_PROJECT="/opt/Emma"

# Copy Emma source code
RUN git clone https://github.com/ian-small/Emma /opt/Emma
        
# Set up Julia environment for Emma
WORKDIR /opt/Emma

# Install Emma dependencies and precompile
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Add additional required packages to the Emma project
RUN julia --project=. -e 'using Pkg; \
        Pkg.add("FASTX"); \
        Pkg.add("BioSequences")'

# Precompile Julia packages and create cache files
RUN julia --project=. -e 'using Pkg; Pkg.precompile()'

# Load Emma to create all necessary cache files
RUN julia --project=. -e 'using Emma; println("Emma loaded successfully")'

# Copy extract_proteins script
COPY extract_proteins.jl /opt/.
RUN chmod +x /opt/*.jl

# Set working directory back to root
WORKDIR /

# Final test
RUN julia -e 'using Emma; println("Emma available in global environment")'




