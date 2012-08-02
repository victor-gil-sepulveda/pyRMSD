#!/bin/bash
cd src

cd cuda
echo "Compiling CUDA RMSD lib"
nvcc -arch=sm_21 --compiler-options '-fPIC' -c RMSDCUDA.cu kernel_rmsd.cu -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart
cd ..

cd serial
echo "Compiling C extensions"
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDSerial.cpp 
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSD.cpp 
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDTools.cpp  
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDomp.cpp  -fopenmp
cd ..

cd python
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC -I/usr/lib64/python2.7/site-packages/numpy/core/include -I/usr/include/python2.7 -c NumpyHelperFuncs.cpp  
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC -I/usr/lib64/python2.7/site-packages/numpy/core/include -I/usr/include/python2.7 -c pyRMSD.cpp  
cd ..

gcc -pthread -shared cuda/RMSDCUDA.o cuda/kernel_rmsd.o python/NumpyHelperFuncs.o python/pyRMSD.o serial/RMSDomp.o serial/RMSD.o serial/RMSDSerial.o serial/RMSDTools.o -L/usr/local/cuda/lib64 -L/usr/lib64 -llapack -lgomp -lcudart -lpython2.7 -o pyRMSD_cfuncs.so -fopenmp

cd ..

mv src/pyRMSD_cfuncs.so pyRMSD/

