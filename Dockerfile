
FROM ubuntu:22.04

ENV SRC /usr/local/src/
ENV BIN /usr/local/bin/
ENV HOME /root/

# install app dependencies
RUN apt-get update && apt-get install -y python3.10 python3-pip
RUN apt-get install -y git
RUN apt-get install -y vim
RUN apt-get install -y htop
RUN apt-get install -y wget


WORKDIR $BIN
RUN wget https://github.com/ddelgadillod/ProtegePD/raw/main/muscle/muscle_lin
RUN chmod +x muscle_lin

WORKDIR $SRC
RUN git clone https://github.com/ddelgadillod/ProtegePD
RUN pip3 install -r ProtegePD/requirements.txt
RUN cp ProtegePD/*.py $BIN
RUN cp -r ProtegePD/assets $BIN
RUN cp -r ProtegePD/test_files $HOME

WORKDIR $BIN
RUN ln -s protege.py protege-pd
WORKDIR $HOME







