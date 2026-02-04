# ============================================================
#  pLIN — Plasmid Life Identification Number System
#  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
#
#  Docker container for self-hosted deployment.
#
#  Build:  docker build -t plin .
#  Run:    docker run -p 8501:8501 plin
#  Open:   http://localhost:8501
# ============================================================

FROM python:3.12-slim

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Working directory
WORKDIR /app

# Install Python dependencies first (cached layer)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY plin_app.py .
COPY train_nt_classifier.py .
COPY .streamlit/ .streamlit/

# Copy data directory (classifier models)
COPY data/ data/

# Copy README for the app's about section
COPY README.md .
COPY LICENSE .
COPY CITATION.cff .

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Run Streamlit
ENTRYPOINT ["streamlit", "run", "plin_app.py", \
            "--server.port=8501", \
            "--server.address=0.0.0.0", \
            "--server.headless=true", \
            "--browser.gatherUsageStats=false"]
