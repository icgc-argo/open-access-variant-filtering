FROM ubuntu:20.04

LABEL \
  version="1.11" \
  description="bcftools image for use in Workflows" \
  org.opencontainers.image.source="https://github.com/icgc-argo/variant-calling-tools"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y python3-pip python3-dev software-properties-common curl \
    && apt-get install -y bzip2 g++ make ncurses-dev zlib1g-dev libz-dev pkg-config libbz2-dev liblzma-dev 

RUN mkdir /tools

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.11

WORKDIR /tmp
RUN curl -sSL -o bcftools-$BCFTOOLS_VERSION.tar.bz2 --retry 10 https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 \
  && tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN ./configure --prefix=$BCFTOOLS_INSTALL_DIR \
  && make \
  && make install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools \
  && rm -rf /tmp/bcftools-$BCFTOOLS_VERSION

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

RUN groupadd -g 1000 ubuntu &&\
    useradd -l -u 1000 -g ubuntu ubuntu &&\
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash"]
