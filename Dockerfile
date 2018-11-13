#using neurodebian runtime as parent image
FROM neurodebian:xenial-non-free
MAINTAINER C-PAC team <cnl@childmind.org>

RUN mkdir -p /code 

#Run the only command 
RUN echo Please use the fcp-indi container instead. Use docker pull fcpindi/c-pac!


#entrypoint
ENTRYPOINT = ["/code/run.py"]
