#!/bin/csh -fe

set dir=$1
cd ./batchfiles/batch10/ddfiles/run$dir 

echo "Executing run $1 on" `hostname` "in $PWD"
setenv PATH "/home1/clage/Software/yt-x86_64/bin:/nasa/sles11/gcc/4.4.4/bin:/nasa/gsl/1.14/bin:/nasa/intel/Compiler/11.1/046/bin/intel64:/nasa/intel/Compiler/11.1/046/mkl/tools/environment:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/PBS/bin:/usr/sbin:/sbin:/opt/c3/bin:/opt/sgi/sbin:/opt/sgi/bin:/nasa/hdf5/1.6.5/serial/bin:/nasa/glib/2.22.4/bin"
setenv LD_LIBRARY_PATH "/nasa/glib/2.22.4/lib:/nasa/sles11/gcc/4.4.4/lib:/nasa/sles11/gcc/4.4.4/lib64:/nasa/gsl/1.14/lib:/nasa/hdf5/1.6.5/serial/lib:/nasa/intel/Compiler/11.1/046/lib/intel64:/nasa/intel/Compiler/11.1/046/mkl/lib/em64t:/nasa/intel/Compiler/11.1/046/ipp/em64t/sharedlib:/home1/clage/Software/yt-x86_64/lib/"
#echo $LD_LIBRARY_PATH
setenv PYTHONPATH "/home1/clage/Software/epd-7.2-2-rh5-x86_64/lib/python2.7/site-packages/:/home1/clage/Software/yt-x86_64/src/yt-hg"

/home1/clage/Research/bullet/code/TriaxialHalo_SGas8GTemp/TriaxialHalo triaxialm.cfg >&sgasm.out

/home1/clage/Research/bullet/code/TriaxialHalo_SGas8GTemp/TriaxialHalo triaxialb.cfg >&sgasb.out

/home1/clage/Research/bullet/code/CombineGalaxies_NewIP/CombineGalaxies bullet.dat 0.0 0.0 0.0 main.dat 173.1 47.4 179.1 196.7237 24.3393 2800.0 0.9134 collision.dat >&collision.out

python /home1/clage/Research/bullet/code/pysubs/mod_grid_locations_big.py  

python /home1/clage/Research/bullet/code/pysubs/make_cool_tfudge.py 0.9 0.4 >&cool.out

python /home1/clage/Research/bullet/code/pysubs/translate_gadget_to_enzo_DEF_big.py 6.95 RandomRadial 4.0 -8.0 >&tran.out

