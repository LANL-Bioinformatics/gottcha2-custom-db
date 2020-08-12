# Build this image:  docker build -t mpi .
#

FROM ubuntu:18.04


MAINTAINER Mark Flynn <mflynn@lanl.gov>



RUN apt-get update -y && \
    apt-get install -y --no-install-recommends sudo apt-utils && \
    apt-get install -y --no-install-recommends openssh-server build-essential wget \
        python-dev python-numpy python-pip python-virtualenv python-scipy \
        gcc gfortran libopenmpi-dev openmpi-bin openmpi-common openmpi-doc binutils zlib1g-dev cmake && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
RUN bash Anaconda3-2020.02-Linux-x86_64.sh -b
RUN rm Anaconda3-2020.02-Linux-x86_64.sh
ENV PATH "/scripts:/root/anaconda3/bin:$PATH"
# Updating Anaconda packages
RUN conda init bash &&  conda update conda && conda update anaconda && conda update --all
RUN conda install -yc \
    bioconda pandas requests minimap2 \
    && conda clean -afy

COPY ./scripts /scripts

COPY ./dbgen /dbgen
WORKDIR /dbgen
RUN make

ENV USER mpirun

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/home/${USER}
RUN mkdir /var/run/sshd
RUN echo 'root:${USER}' | chpasswd
RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# ------------------------------------------------------------
# Add an 'mpirun' user
# ------------------------------------------------------------

RUN adduser --disabled-password --gecos "" ${USER} && \
    echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# ------------------------------------------------------------
# Set-Up SSH with our Github deploy key
# ------------------------------------------------------------

ENV SSHDIR ${HOME}/.ssh/

RUN echo ${USER}
RUN echo ${HOME}
RUN mkdir -p ${SSHDIR}

RUN echo ${SSHDIR}
ADD ssh/config ${SSHDIR}/config
ADD ssh/id_rsa.mpi ${SSHDIR}/id_rsa
ADD ssh/id_rsa.mpi.pub ${SSHDIR}/id_rsa.pub
ADD ssh/id_rsa.mpi.pub ${SSHDIR}/authorized_keys

RUN chmod -R 600 ${SSHDIR}* && \
    chown -R ${USER}:${USER} ${SSHDIR}

RUN pip install --upgrade pip

USER ${USER}
#RUN  pip install --user -U setuptools \
#    && pip install --user mpi4py

# ------------------------------------------------------------
# Configure OpenMPI
# ------------------------------------------------------------

USER root

RUN rm -fr ${HOME}/.openmpi && mkdir -p ${HOME}/.openmpi
ADD default-mca-params.conf ${HOME}/.openmpi/mca-params.conf
RUN chown -R ${USER}:${USER} ${HOME}/.openmpi

# ------------------------------------------------------------
# Copy MPI4PY example scripts
# ------------------------------------------------------------

#ENV TRIGGER 1

#ADD mpi4py_benchmarks ${HOME}/mpi4py_benchmarks
#RUN chown -R ${USER}:${USER} ${HOME}/mpi4py_benchmarks


EXPOSE 22
#CMD ["/usr/sbin/sshd", "-D"]
ENTRYPOINT ["/bin/bash","-c"]
CMD ["/bin/bash"]