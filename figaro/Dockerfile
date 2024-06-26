FROM python:3.6.7

MAINTAINER Michael M. Weinstein, Zymo Research
LABEL version="1.0.0"

WORKDIR /opt

# Clonar el repositorio
RUN git clone https://github.com/Zymo-Research/figaro.git

# Navegar al directorio del proyecto
WORKDIR /opt/figaro

# Instalar las dependencias
RUN pip3 install --upgrade pip && \
    pip3 install -r requirements.txt

ENV PYTHONUNBUFFERED=1

CMD ["python3", "/opt/figaro/figaro/figaro.py", "-h"]
