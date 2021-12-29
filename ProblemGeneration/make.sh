#!/bin/bash
#

#download the files
mkdir dataSet
cd ./dataSet

mkdir real-sim
cd real-sim
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/real-sim.bz2
bunzip2 real-sim.bz2
mkdir SVMtoLP
cd ..

mkdir news20
cd news20
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass/news20.bz2
bunzip2 news20.bz2
mkdir SVMtoLP
cd ..

mkdir rcv1
cd rcv1
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass/rcv1_train.multiclass.bz2
bunzip2 rcv1_train.multiclass.bz2
mv rcv1_train.multiclass rcv1
mkdir SVMtoLP
cd ..

mkdir avazu
cd avazu
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/avazu-app.bz2
bunzip2 avazu-app.bz2
mv avazu-app avazu
mkdir SVMtoLP
cd ..

#compile
cd ./Problem
g++ -std=c++11 -o SVMtoLP SVMtoLP.cpp

#generate SVMtoLP problems
./SVMtoLP -t 1 -k 2 real-sim
./SVMtoLP -t 1 -k 20 news20
./SVMtoLP -t 1 -k 51 rcv1
./SVMtoLP -t 1 -k 2 avazu
