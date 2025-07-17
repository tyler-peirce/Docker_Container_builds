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

# Create directories
RUN mkdir -p /opt/julia-depot /tmp/julia-depot /tmp/julia/pkgcache /tmp/julia/load_cache /tmp/julia

# Add Julia to PATH
ENV PATH="/opt/julia-1.10.5/bin:$PATH"
# Set Julia load path to include Emma
ENV JULIA_LOAD_PATH="/opt/Emma:/opt/Emma/src:$JULIA_LOAD_PATH"
# Make sure the Julia environment is activated by default
ENV JULIA_PROJECT="/opt/Emma"
ENV JULIA_DEPOT_PATH="/tmp/julia-depot"
ENV JULIA_HISTORY="/tmp/julia/repl_history.jl"
ENV JULIA_PKG_CACHE="/tmp/julia/pkgcache"
ENV JULIA_LOAD_CACHE_PATH="/tmp/julia/load_cache"
ENV HOME="/opt"
ENV JULIA_PKG_PRECOMPILE_AUTO=0


# Copy Emma source code
RUN git clone https://github.com/ian-small/Emma /opt/Emma
        
# Set up Julia environment for Emma
WORKDIR /opt/Emma

# Add additional required packages to the Emma project
RUN julia --project=. -e 'using Pkg; \
        Pkg.instantiate(); \
        Pkg.add("FASTX"); \
        Pkg.add("BioSequences"); \
        Pkg.precompile(); \
        using Emma'

# Load Emma to create all necessary cache files
RUN julia --project=. -e 'using Emma; println("Emma loaded successfully")'

# Copy extract_proteins script
COPY extract_proteins.jl /opt/.
RUN chmod +x /opt/*.jl
RUN mv /opt/julia-depot /tmp/julia-depot
# Set working directory back to root
WORKDIR /

# Final test
RUN julia -e 'using Emma; println("Emma available in global environment")'




