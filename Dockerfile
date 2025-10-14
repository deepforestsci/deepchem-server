FROM mambaorg/micromamba:2.3.2

ARG DEEPCHEM_SERVER_HOME=/opt/deepchem_server_app
ARG MAMBA_USER=mambauser

ENV DEEPCHEM_SERVER_HOME=${DEEPCHEM_SERVER_HOME} \
    MAMBA_NO_LOW_SPEED_LIMIT=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=random \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    DEBIAN_FRONTEND=noninteractive \
    DATADIR=${DEEPCHEM_SERVER_HOME}/data

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    libgomp1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
USER $MAMBA_USER

WORKDIR ${DEEPCHEM_SERVER_HOME}

COPY --chown=$MAMBA_USER:$MAMBA_USER deepchem_server/environments/core_environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes && \
    rm -f /tmp/env.yaml

RUN mkdir -p ${DATADIR} && \
    chmod 755 ${DATADIR}

COPY --chown=$MAMBA_USER:$MAMBA_USER deepchem_server/ ${DEEPCHEM_SERVER_HOME}/deepchem_server/

WORKDIR ${DEEPCHEM_SERVER_HOME}

EXPOSE 8000

HEALTHCHECK --interval=120s --timeout=10s --start-period=180s --retries=3 \
    CMD curl -f http://127.0.0.1:8000/healthcheck || exit 1

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]