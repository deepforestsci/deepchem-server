FROM continuumio/miniconda3:23.5.2-0

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app/deepchem_server

COPY deepchem_server/environments/core_environment.yml ./

# Install system dependencies and clean up in a single layer
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    curl \
    libgomp1 \
    git \
    libboost-all-dev \
    swig \
    build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    conda update -n base -c defaults conda && \
    conda env update -n base -f ./core_environment.yml && \
    conda clean -afy

COPY deepchem_server/ .

EXPOSE 8000

WORKDIR /app

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/healthcheck || exit 1

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
