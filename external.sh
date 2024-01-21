#!/bin/bashi
echo 'Download muscle'
cd ProtegePD
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz
tar -xzvf muscle3.8.31_i86linux32.tar.gz
rm muscle3.8.31_i86linux32.tar.gz
chmod +x muscle3.8.31_i86linux32
mkdir -p ProtegePD/bin
mv muscle3.8.31_i86linux32 ./bin/muscle
./bin/muscle
cd ..

