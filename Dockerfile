FROM python:3.13.1-slim-bookworm

WORKDIR /opt

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY scripts/* .
COPY resources/* .

CMD [ "python", "/opt/Parse-MCI_JSONs.py" ]
