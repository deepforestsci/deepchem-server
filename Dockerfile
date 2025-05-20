FROM python:3.10-slim

WORKDIR /app/deepchem_server

COPY deepchem_server/requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY deepchem_server/ .

EXPOSE 8000

WORKDIR /app

CMD ["uvicorn", "deepchem_server.main:app", "--host", "0.0.0.0", "--port", "8000"]
