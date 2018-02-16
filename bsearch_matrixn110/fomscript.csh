#!/bin/csh -f

@ rank = 1 + $MPIEXEC_RANK

cd ./batchfiles/batch10/ddfiles/run$rank 
python /home1/clage/Research/bullet/code/pysubs/findfom_26Sep13.py >&fom.out
