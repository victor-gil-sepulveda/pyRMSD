#!/bin/bash
cd src

cd theobald
echo -e "\033[31mCompiling CUDA RMSD lib\033[0m"

nvcc --use_fast_math -arch sm_11 --compiler-options '-fPIC' -I/usr/local/cuda-4.2/include -c ThRMSDCuda.cu kernel_functions_cuda.cu  

gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c ThRMSDSerial.cpp

gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c kernel_functions_serial.cpp

gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c ThRMSDSerialOmp.cpp

gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c kernel_functions_omp.cpp -fopenmp

	cd test/
	echo -e "\033[31m\tCompiling C test\033[0m"
	nvcc -I/usr/local/cuda-4.2/include -c ctest.cpp 
	gcc -c ../../serial/RMSD.cpp 
	g++ ctest.o RMSD.o  ../ThRMSDSerial.o ../kernel_functions_serial.o ../ThRMSDCuda.o ../kernel_functions_cuda.o  -L/usr/local/cuda-4.2/lib64 -octest -lgomp -lcudart
	cd ..
cd ..

cd serial
echo -e "\033[31mCompiling Serial and OpenMP C lib\033[0m"
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDSerial.cpp 
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSD.cpp 
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDTools.cpp  
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c RMSDomp.cpp  -fopenmp
cd ..

cd matrix
echo -e "\033[31mCompiling Matrix lib\033[0m"
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O3 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables  -fPIC  -c matrix.cpp 
cd ..

cd python
echo -e "\033[31mCompiling Python bindings\033[0m"
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/python2.7 -c NumpyHelperFuncs.cpp  
gcc -pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/python2.7 -I/usr/local/cuda-4.2/include -c pyRMSD.cpp  
cd ..

echo -e "\033[31mLinking...\033[0m"
g++ -pthread -shared  theobald/ThRMSDSerial.o theobald/ThRMSDSerialOmp.o theobald/ThRMSDCuda.o theobald/kernel_functions_serial.o theobald/kernel_functions_omp.o theobald/kernel_functions_cuda.o python/NumpyHelperFuncs.o python/pyRMSD.o serial/RMSDomp.o serial/RMSD.o serial/RMSDSerial.o serial/RMSDTools.o matrix/matrix.o -L/usr/local/cuda-4.2/lib64 -llapack -lgomp -lcudart -lpython2.7 -o pyRMSD_cfuncs.so -fopenmp

g++ -pthread -shared matrix/matrix.o -lpython2.7 -o condensedMatrix.so -fopenmp
cd ..

echo -e "\033[31mMoving...\033[0m"
mv src/pyRMSD_cfuncs.so pyRMSD/
mv src/condensedMatrix.so pyRMSD/
