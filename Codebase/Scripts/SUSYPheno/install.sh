# This script compiles all the needed packages
export SUSYPHENO_PATH=$PWD
echo $SUSYPHENO_PATH
cd $SUSYPHENO_PATH/HDECAY
make clean
make
cd $SUSYPHENO_PATH/MICROMEGA/micromegas_3.5.5/
make clean
mv CalcHEP_src/FlagsForSh CalcHEP_src/FlagsForSh_obso
gmake
cd $SUSYPHENO_PATH/MICROMEGA/micromegas_3.5.5/MSSM/
gmake main=micromegas_mini_slha.c
cd $SUSYPHENO_PATH/SUSYHIT/
make clean
make
cd $SUSYPHENO_PATH/DARKSUSY/darksusy-5.1.1/
make uninstall
make distclean
./configure
make
make install
cd $SUSYPHENO_PATH
mkdir $SUSYPHENO_PATH/MICROMEGA/micromegas_3.5.5/MSSM/work/so_generated/
export PATH=$PATH:$SUSYPHENO_PATH/bin

