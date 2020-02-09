FROM ubuntu
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update -y
RUN apt-get install -y libstdc++6 wget build-essential libx11-dev libboost-dev libcurl4 cmake swig flex bison libtool zlib1g-dev tcl-dev tk-dev libjpeg-dev cpio
# RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
# RUN sh -c 'echo deb https://apt.repos.intel.com/ipp all main > /etc/apt/sources.list.d/intel-ipp.list'
# RUN apt-get update -y
# RUN apt-get install -y intel-mkl-2018.2-046 intel-ipp-2018.4-057 
RUN wget http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16318/l_mkl_2020.0.166.tgz
RUN tar -xvf l_mkl_2020.0.166.tgz
RUN rm l_mkl_2020.0.166.tgz
RUN l_mkl_2020.0.166/install.sh -s l_mkl_2020.0.166/silent.cfg --accept_eula
RUN rm -r l_mkl_2020.0.166
RUN wget http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16233/l_ipp_2020.0.166.tgz
RUN tar -xvf l_ipp_2020.0.166.tgz
RUN rm l_ipp_2020.0.166.tgz
RUN l_ipp_2020.0.166/install.sh -s l_ipp_2020.0.166/silent.cfg --accept_eula
RUN rm -r l_ipp_2020.0.166
RUN mkdir /preplace
RUN apt-get install -y unzip exuberant-ctags
ADD bookshelf.zip /
RUN unzip bookshelf.zip
ENV PREPLACEPATH /preplace
ENTRYPOINT [ "/preplace/src/bookshelf.sh", "1", "adaptec1", "f.csv" ]
ADD . /preplace/
RUN mv bookshelf /preplace/bench/bookshelf
#RUN make hash -C /preplace
#RUN make prep -C /preplace
RUN make -C /preplace