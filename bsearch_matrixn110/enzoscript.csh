#!/bin/csh -fe

set cpus_per_node=8
set nodes_per_job=8
@ cpus_per_job=$cpus_per_node * $nodes_per_job

set dir=$1

if ($dir % $nodes_per_job == 0) then
    @ dir=$dir / $nodes_per_job
    cd ./batchfiles/batch10/ddfiles/run$dir 
    set count=0
    sleep 20.0
    while ($count < $cpus_per_node)
	echo `hostname` >> ../hostfile_$dir
	@ count += 1
    end
    sleep 20.0
    echo "Executing run $1 on" `hostname` "in $PWD"
    setenv PBS_NODEFILE ../hostfile_$dir
    set path = (/nasa/sgi/mpt/2.04.10789/bin $PATH)
    setenv LD_LIBRARY_PATH "/nasa/sgi/mpt/2.04.10789/lib:/nasa/glib/2.22.4/lib:/nasa/sles11/gcc/4.4.4/lib:/nasa/sles11/gcc/4.4.4/lib64:/nasa/gsl/1.14/lib:/nasa/hdf5/1.6.5/serial/lib:/nasa/intel/Compiler/11.1/046/lib/intel64:/nasa/intel/Compiler/11.1/046/mkl/lib/em64t:/nasa/intel/Compiler/11.1/046/ipp/em64t/sharedlib"
    set first_cpu=0
    set last_cpu=$cpus_per_node
    @ last_cpu -= 1
    mpiexec -np $cpus_per_job dplace -s1 -c$first_cpu-$last_cpu /home1/clage/Software/enzo-2.0/src/enzo_MHD_Visc14/enzo.exe AMRTest.enzo >&enzo.out&

else if ($dir % $nodes_per_job == 1) then
    @ dir=$dir - $dir % $nodes_per_job 
    @ dir=$dir / $nodes_per_job + 1
    cd ./batchfiles/batch10/ddfiles/run$dir 
    echo `hostname` >! ../hostfile_$dir
    set count=1
    sleep 20.0
    while ($count < $cpus_per_node)
	echo `hostname` >> ../hostfile_$dir
	@ count += 1
    end

else
    @ dir=$dir - $dir % $nodes_per_job 
    @ dir=$dir / $nodes_per_job + 1
    cd ./batchfiles/batch10/ddfiles/run$dir 
    set count=0
    sleep 20.0
    while ($count < $cpus_per_node)
	echo `hostname` >> ../hostfile_$dir
	@ count += 1
    end
endif
wait
