FROM python:3.10

WORKDIR /app/deepchem_server

COPY deepchem_server/requirements.txt .

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get -y install curl libgomp1

RUN pip install --no-cache-dir -r requirements.txt

COPY deepchem_server/ .

EXPOSE 8000

WORKDIR /app

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
