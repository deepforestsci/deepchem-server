FROM python:3.10-slim

WORKDIR /app

RUN pip install --no-cache-dir fastapi uvicorn httpx python-multipart

COPY deepchem_server /app/deepchem_server

ENV DATASTORE_BASE_DIR=/data
ENV DATASTORE_PORT=8081
ENV DATASTORE_HOST=0.0.0.0

RUN mkdir -p /data

EXPOSE 8081

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import httpx; httpx.get('http://localhost:8081/api/v1/healthcheck').raise_for_status()"

CMD ["python", "-m", "deepchem_server.services.datastore"]
