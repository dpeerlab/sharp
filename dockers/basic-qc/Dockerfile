FROM sailmskcc/hto-dnd:1.1.2

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

COPY render /opt/render
COPY templates /opt/templates


ENTRYPOINT [ "python", "/opt/render/cli.py" ]