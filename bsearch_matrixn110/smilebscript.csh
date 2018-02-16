#!/bin/csh -fe

set dir=$1
cd ./batchfiles/batch2/ddfiles/run$dir 

echo "Executing run $1 on" `hostname` "in $PWD"

#set path = (/nasa/sgi/mpt/2.04.10789/bin $PATH)
#echo $PATH
setenv LD_LIBRARY_PATH "/nasa/glib/2.22.4/lib:/nasa/sles11/gcc/4.4.4/lib:/nasa/sles11/gcc/4.4.4/lib64:/nasa/gsl/1.14/lib:/nasa/hdf5/1.6.5/serial/lib:/nasa/intel/Compiler/11.1/046/lib/intel64:/nasa/intel/Compiler/11.1/046/mkl/lib/em64t:/nasa/intel/Compiler/11.1/046/ipp/em64t/sharedlib"
#echo $LD_LIBRARY_PATH
/home1/clage/Research/bullet/code/TriaxialHalo_Smile203/smilec smileb.script >&smileb.out

