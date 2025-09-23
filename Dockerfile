FROM mambaorg/micromamba:2.3.2

ENV MAMBA_NO_LOW_SPEED_LIMIT=1 \
	PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app/deepchem_server

COPY --chown=$MAMBA_USER:$MAMBA_USER deepchem_server/environments/core_environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

COPY deepchem_server/ .

EXPOSE 8000

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgomp1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

WORKDIR /app

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/healthcheck || exit 1

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
