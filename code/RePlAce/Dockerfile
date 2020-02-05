FROM ubuntu
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update -y
RUN apt-get install -y libstdc++6 wget build-essential libx11-dev libboost-dev libcurl4 cmake swig flex bison libtool zlib1g-dev tcl-dev tk-dev libjpeg-dev
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN sh -c 'echo deb https://apt.repos.intel.com/ipp all main > /etc/apt/sources.list.d/intel-ipp.list'
RUN apt-get update -y
RUN apt-get install -y intel-mkl-2018.2-046 intel-ipp-2018.4-057 
RUN mkdir /preplace
RUN apt-get install -y unzip
RUN wget http://download2263.mediafire.com/4adno2bnm6jg/n6jkth9a3zpd9nm/bookshelf.zip 
RUN unzip bookshelf.zip
RUN apt-get install -y exuberant-ctags
ENV PREPLACEPATH /preplace
ENTRYPOINT [ "/preplace/src/bookshelf.sh", "1", "adaptec1", "f.csv" ]
ADD . /preplace/
RUN mv bookshelf /preplace/bench/bookshelf
#RUN make hash -C /preplace
#RUN make prep -C /preplace
RUN make -C /preplace