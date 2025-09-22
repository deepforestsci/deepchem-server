FROM mambaorg/micromamba:1.4.8

WORKDIR /app/deepchem_server

ENV MAMBA_NO_LOW_SPEED_LIMIT=1

COPY --chown=$MAMBA_USER:$MAMBA_USER deepchem_server/environments/core_environment.yml ./core_environment.yml

USER root
RUN apt update && apt install -y --no-install-recommends apt-utils && apt install -y build-essential git wget curl libgomp1
USER $MAMBA_USER

RUN micromamba install -y -n base -f ./core_environment.yml

COPY deepchem_server/ .

EXPOSE 8000

WORKDIR /app

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
