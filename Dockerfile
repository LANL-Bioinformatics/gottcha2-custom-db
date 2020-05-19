#FROM nlknguyen/alpine-mpich
FROM sourceryinstitute/mpich-docker-image:latest
#RUN sudo apk add zlib-dev
RUN apt-get update && apt-get install -y build-essential zlib1g-dev && \
    apt-get autoremove -y && apt-get clean -y
COPY ./ ./

RUN make
ENTRYPOINT ["/bin/bash","-c"]
CMD ["bin/bash -l"]