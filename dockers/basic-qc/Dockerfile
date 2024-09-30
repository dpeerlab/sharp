FROM python:3.12-slim-bullseye

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

RUN pip install nbconvert==5.6.1 jupyter_contrib_nbextensions==0.5.1 jinja2==3.0.3
RUN jupyter contrib nbextension install --user \
&& jupyter nbextension enable toc2
RUN python -m ipykernel install --name=python_qc

COPY ./notebooks/inspect-*.ipynb /opt/
