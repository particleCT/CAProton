module add cmake/2.8.12.2
source 
export G4WORKDIR=$PWD
#export LD_LIBRARY_PATH=/home/cacof1/geant4.9.6.p03/lib64/:$LD_LIBRARY_PATH
mkdir build
cd build

rm -rf *


cmake -DGeant4_DIR=/Users/cacof1/software/geant4.10.00/lib/Geant4-10.0.3 -DCMAKE_INSTALL_PREFIX=$G4WORKDIR ../

make -j8

make install

cd ../
rm -rf ./build/
#./bin/Darwin-clang/CAProton 100000 $RANDOM Slab 300 10