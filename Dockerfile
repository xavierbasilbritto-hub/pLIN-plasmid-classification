# ============================================================
#  pLIN — Plasmid Lineage Identification Number System
#  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
#
#  Docker container with all bioinformatics tools pre-installed.
#  Works on Linux, macOS (Intel/Apple Silicon), and Windows (Docker Desktop).
#
#  Build:  docker build -t plin .
#  Run:    docker run -p 8501:8501 plin
#  Open:   http://localhost:8501
#
#  Full build with bioinformatics tools (recommended for reviewers):
#  docker build --build-arg INSTALL_BIOTOOLS=true -t plin-full .
# ============================================================

FROM python:3.11-slim AS base

LABEL maintainer="Basil Xavier Britto"
LABEL description="pLIN: Plasmid Lineage Identification Number System"
LABEL version="2.1.0"

# Build argument: set to "true" to install bioinformatics tools
ARG INSTALL_BIOTOOLS=true

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# ── System dependencies ─────────────────────────────────────────────────
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        wget \
        ca-certificates \
        libgomp1 \
        procps \
        default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# ── Install Miniconda (for bioconda tools) ───────────────────────────────
RUN if [ "$INSTALL_BIOTOOLS" = "true" ]; then \
        wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
        bash /tmp/miniconda.sh -b -p /opt/conda && \
        rm /tmp/miniconda.sh && \
        /opt/conda/bin/conda config --add channels defaults && \
        /opt/conda/bin/conda config --add channels bioconda && \
        /opt/conda/bin/conda config --add channels conda-forge && \
        /opt/conda/bin/conda config --set channel_priority strict && \
        /opt/conda/bin/conda clean -afy; \
    fi

ENV PATH="/opt/conda/bin:${PATH}"

# ── Install bioinformatics tools via conda ───────────────────────────────
RUN if [ "$INSTALL_BIOTOOLS" = "true" ]; then \
        conda install -y -c bioconda -c conda-forge \
            ncbi-amrfinderplus \
            mash \
            fastani \
            minimap2 \
            minced \
            blast \
            prodigal \
        && conda clean -afy \
        && amrfinder --update 2>/dev/null || true; \
    fi

# ── Working directory ────────────────────────────────────────────────────
WORKDIR /app

# ── Install Python dependencies (cached layer) ──────────────────────────
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ── Copy application code ────────────────────────────────────────────────
COPY plin_app.py .
COPY assign_pLIN.py .
COPY assign_pLIN_reference.py .
COPY build_inc_centroids.py .
COPY integrate_pLIN_AMR.py .

# ── Copy data directory (classifier models, ~5 MB) ──────────────────────
COPY data/ data/

# ── Copy output directory (pLIN assignments TSV for query mode) ──────────
COPY output/pLIN_assignments.tsv output/pLIN_assignments.tsv

# ── Copy Streamlit config ────────────────────────────────────────────────
COPY .streamlit/ .streamlit/

# ── Copy documentation ───────────────────────────────────────────────────
COPY README.md .
COPY LICENSE .
COPY CITATION.cff .

# ── Create output directories ────────────────────────────────────────────
RUN mkdir -p output/amrfinder output/crispr_evaluation

# ── Expose Streamlit port ────────────────────────────────────────────────
EXPOSE 8501

# ── Health check ─────────────────────────────────────────────────────────
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# ── Run Streamlit ────────────────────────────────────────────────────────
ENTRYPOINT ["streamlit", "run", "plin_app.py", \
            "--server.port=8501", \
            "--server.address=0.0.0.0", \
            "--server.headless=true", \
            "--server.maxUploadSize=200", \
            "--browser.gatherUsageStats=false"]
