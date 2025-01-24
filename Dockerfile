FROM ubuntu:focal
RUN apt update && apt install -y python3 python3-pip python3-gmsh
RUN pip3 install numpy
RUN pip3 install matplotlib
RUN pip3 install pygmsh
