FROM ubuntu:22.04
RUN apt update -y
RUN apt install -y gnupg ca-certificates
RUN apt-key adv --fetch-key https://repo.arrayfire.com/GPG-PUB-KEY-ARRAYFIRE-2020.PUB
RUN echo "deb [arch=amd64] https://repo.arrayfire.com/debian all main" | tee /etc/apt/sources.list.d/arrayfire.list
RUN apt update -y
RUN DEBIAN_FRONTEND=noninteractive apt install -y build-essential gnuradio libiio* cmake arrayfire-unified3-dev intel-mkl arrayfire-opencl3-dev arrayfire-opencl3-mkl arrayfire-unified3 libbtbb-dev pocl-opencl-icd
RUN DEBIAN_FRONTEND=noninteractive apt install -y git
RUN git clone https://github.com/CellWizard/gr-bluetooth-plutoblue.git
RUN git clone https://github.com/bkerler/gr-bluetooth.git
RUN cd gr-bluetooth && cmake . && make install && cd ..
RUN apt install -y arrayfire-cpu3-mkl vim
ADD apps /apps
ADD system_top.bit.bin /system_top.bit.bin
ADD run.sh /run.sh
add loadplutobluez_runme_on_pluto.txt /loadplutobluez_runme_on_pluto.txt
ADD plutoblueconverter /plutoblueconverter

