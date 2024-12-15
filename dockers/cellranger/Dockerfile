FROM almalinux/amd64:9.3-base-20231124

# CentOS 8 is dead, so we are using AlmaLinux 8.9. Built by Oracle.

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

ARG CELLRANGER_VERSION
ARG DOWNLOAD_URL
ENV PATH /opt/cellranger-${CELLRANGER_VERSION}:$PATH

RUN yum group install -y "Development Tools" \
    && yum install -y which

# cell ranger binaries
RUN curl -o cellranger-${CELLRANGER_VERSION}.tar.gz ${DOWNLOAD_URL} \
    && tar xzf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && rm -rf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && mv cellranger-${CELLRANGER_VERSION} /opt/

# Removing VDJ reference data compared to previous version to make the container smaller

WORKDIR /opt

ENTRYPOINT [ "cellranger" ]
