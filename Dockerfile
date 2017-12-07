FROM python:2.7-alpine

# Download and install samtools:
RUN apk add --no-cache --virtual build-dependencies alpine-sdk \
	&& curl -L https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 -o samtools.tar.bz2 \
	&& tar -jxf samtools.tar.bz2 \
	&& cd samtools-1.6 \
	&& mkdir /usr/local/samtools \
	&& ./configure --prefix=/usr/local/samtools \
	&& make \
	&& make -install \
	&& cd \
	&& ls -l \
	&& apk del build-dependencies