
FROM ubuntu:22.04

ENV SRC /usr/local/src/
ENV BIN /usr/local/bin/

# install app dependencies
RUN apt-get update && apt-get install -y python3.10 python3-pip
RUN apt-get install -y git
RUN apt-get install -y vim
RUN apt-get install -y htop

#WORKDIR $BIN
#RUN wget https://github.com/ddelgadillod/ProtegePD/raw/main/muscle/muscle_lin

#WORKDIR $SRC
RUN git clone https://github.com/ddelgadillod/ProtegePD
#RUN pip3 install -r ProtegePD/requirements.txt
#RUN cp *.py $BIN
#RUN cp -r assets $BIN



#RUN cd ProtegePD
#RUN cp muscle/muscle_lin $BIN


# CMD bash install.sh


