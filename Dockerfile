FROM python:2.7

# Download and install samtools:
RUN apt-get update \
	&& apt-get install -y --no-install-recommends build-essential  \
	&& curl -L https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 -o samtools.tar.bz2 \
	&& tar -jxf samtools.tar.bz2 \
	&& cd samtools-1.6 \
	&& mkdir /usr/local/samtools \
	&& ./configure --prefix=/usr/local/samtools \
	&& make \
	&& make -install \
	&& cd \
	&& ls -l \
	
RUN pip install --upgrade pip \
  && pip install -r requirements.txt
  
ENTRYPOINT ["python"]
CMD ["app.py"]
	
	