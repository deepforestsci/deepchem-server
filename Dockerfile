FROM continuumio/miniconda3:23.5.2-0

WORKDIR /app/deepchem_server

COPY deepchem_server/requirements.txt .

COPY deepchem_server/environments/core_environment.yml ./core_environment.yml

RUN apt-get update  \
    && apt-get install -y --no-install-recommends apt-utils  \
    && apt-get install -y curl libgomp1 git libboost-all-dev swig build-essential

# Update conda and install Python 3.11, openmm and pdbfixer from conda-forge
RUN conda update -n base -c defaults conda \
    && conda env update -n base -f ./core_environment.yml \
    && conda clean -ay

COPY deepchem_server/ .

EXPOSE 8000

WORKDIR /app

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
