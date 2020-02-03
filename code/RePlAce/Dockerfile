FROM ubuntu
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update -y
RUN apt-get install -y libstdc++6 wget build-essential libx11-dev libboost-dev libcurl4 cmake swig flex bison libtool zlib1g-dev tcl-dev tk-dev libjpeg-dev
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN sh -c 'echo deb https://apt.repos.intel.com/ipp all main > /etc/apt/sources.list.d/intel-ipp.list'
RUN apt-get update -y
#RUN apt-get install intel-mkl-64bit-2018.2-046
RUN apt-get install -y intel-mkl-2018.2-046 intel-ipp-2018.4-057
RUN mkdir /preplace
ADD . /preplace/
RUN make -C /preplace