# pyRMSD
pyRMSD goal is the fast (and easy!) calculation of rmsd collective operations, specially matrices of large ensembles of protein conformations. It also offers a symmetric distance matrix implementation with improved access speed and memory efficiency.

# Index
>##[1- Features](#features)
>##[2- Building & Installation](#building)  
>>### [Dependencies](#dependencies)  
>>### [Linux](#ilinux)  
>>### [Windows](#iwindows)  
>>### [MacOs](#imac)
>##[3- The custom building script](#buildscript)
>>### [Unix-Based](#build_linux)  
>>### [Windows](#build_win)  
>##[4- Testing (Developers)](#testing)
>##[5- Benchmarcks (Developers)](#benchmarks) 
>##[6- Using it](#usage)
 
##<a id="features"></a> 1- Features
##<a id="building"></a> 2- Building & Installation
### <a id="dependencies"></a> Dependencies
**Users** only need to install Python version 2.6/2.7 (pyRMSD has only been tested with those, however it may work with another versions of the Python 2.X family). Numpy is also required. Surely you have already installed it, but in the case you didn't it can be found [here](http://sourceforge.net/projects/numpy/files/) where you will be able to find installers for almost all the combinations of platforms and Python versions you can think about.

As a **Developer** you may be interested on istalling [scipy](http://www.scipy.org/) (only necessary to execute the statistics test), and [prody](http://www.csb.pitt.edu/prody/getprody.html) (to fully satisfy the RMSD calculators test). We provide our own pdb reader, but only for the sake of completeness. That's why we encourage the use of Prody to handle coordinates, as it is well-tested and powerful tool. Also, remember that header files of Python and Numpy may be accessible, and your Python installation must contain the python shared library This usually means to use ./configure --enable-shared before building Python (usually 2.7 distributions already come with this library).

### <a id="ilinux"></a>Linux 
Linux users have the following choices:  
  
**1)** Using the 'setup.py' file inside the root folder by typing:  
    \> python setup.py install
  
(or 'build') to only build it), which is the usual way python packages are deployed. AS 'distutils' do not support CUDA compiling directly, your package will not be able to use CUDA calculators.

**2)** Using [pyRMSD-1.0.tar.gz](https://github.com/victor-gil-sepulveda/pyRMSD/tree/master/prebuilt_packages/v1/Linux/64-CUDA) precompiled distribution: 

- Unzip the file with:
    \> tar -zxvf
  
- Open pyRMSD folder and use setup.py:
    \> python setup.py install

This distribution will only work in x64 systems, and will make the CUDA calculator available.

**3)** Using the custom build.py script in pyRMSD main folder with:  
    \> python build.py  

or  
    \> python build.py --cuda  

The build.py script is the most versatile way and will work in almost all situations, but as it requires the 'sysconfig' package, the script itself needs Python 2.7 to be executed. With this script one can build x86 and x64 distributions with enabled CUDA calculators if the --cuda flag is used. There is more information about the build.py script [here](#build_linux). 

### <a id="iwindows"></a>Windows
Windows users have the following choices:  
**1)** Using a precompiled windows installer. There are two available ([x86](https://github.com/victor-gil-sepulveda/pyRMSD/tree/master/prebuilt_packages/v1/Win/32) and [x64]()). Those are used as any other regular Windows Installer (double click the executable and follow instructions).  
**2)** Using the custom build.py script in pyRMSD main folder with:  
    \> python build_windows.py  
Please look [here](#build_win) if you need further iformation about the windows version of the custom build script.

### <a id="imac"></a>MacOs
MacOs users have the same choices that Linux users:

**1)** Use the 'setup.py' file inside the root folder by typing:  
    \> python setup.py install
  
(or 'build') to only build it), which is the usual way python packages are deployed. AS 'distutils' do not support CUDA compiling directly, your package will not be able to use CUDA calculators. 

**2)** Using [pyRMSD-1.0.tar.gz](https://github.com/victor-gil-sepulveda/pyRMSD/tree/master/prebuilt_packages/v1/MacOs) precompiled distribution: 

- Unzip the file with:
    \> tar -zxvf
  
- Open pyRMSD folder and use setup.py:
    \> python setup.py install

This is a precompiled distribution without CUDA calculators.

**3)** Using the custom build.py script in pyRMSD main folder with:  
    \> python build.py  

or  
    \> python build.py --cuda  

Please see this same section in the Linux instalation guide. It has only been tested without CUDA support, but it may need only [minor changes](#build_linux) to add it.  

##<a id="buildscript"></a>3- The custom building script  
pyRMSD includes a small build script that is indeed a recipe to compile the C extensions of pyRMSD. This script only works with Python 2.7+ as it uses the module *sysconfig* to get the search path for python headers and libs. 
The building script will try to guess the location of the needed files for compilation, however it can be modified to be able to handle all kind of scenarios. In order to do this you may want to modify the upper case constants in the top part of the file. 

###<a id="build_linux"></a>Unix-based systems  
The script was used in a Ubuntu x86 and Ubuntu x64 Os, as well as a MacOs (Snow Leopard) for the non CUDA build. PYTHON_X constants were left unchanged.
It was also used under Ubuntu x64 with CUDA 4.2 to build the CUDA enabled version. 
If you are going to use it to build a CUDA enabled version you may have to change the *CUDA_BASE* constant, which needs to point to the base directory of your CUDA installation (in our case  */usr/local/cuda-4.2*). Required headers and libs are usually stored inside the */include* and */lib64* folders (*/lib* in x86 systems) subfolders, but you can also change it by modifying *CUDA_INCLUDE_FOLDER* and *CUDA_LIBRARIES_FOLDER*.
Finally you can change *CUDA_ARCHITECHTURE* to match the architecture of your GPU. 
###<a id="build_win"></a>Windows systems  

##<a id="testing"></a>Testing (Developers)  
Once installed you can run the tests in *pyRMSD/test* using:  
  
    \> python -m unittest testCondensedMatrix testMatrixHandler testMatrixNeighbours testMatrixStatistics testRMSDCalculators testPdbReader

Currently only the *test_create_with_reader* test will fail if all the dependencies are fullfilled (it's unwritten yet). 
If you didn't build pyRMSD with CUDA support, 3 tests will be skipped.

##Benchmarcks (Developers)

### Building
To build the code (over Linux, no Windows support yet) you must execute *install.py*.  
  
    > python install.py
  
This will build all serial and OpenMP calculators. *install.py* can be modified to add/remove compiling options (if you know what you're doing).  
To add the CUDA calculators, just use:  
  
    > python install.py --cuda
  
and it will also try to build the available CUDA calculators. In this case you will surely have to modify the *CUDA_BASE*, *CUDA_INCLUDE_FOLDER*, *CUDA_LIB_FOLDER*, *CUDA_ARCHITECHTURE* and *CUDA_LIBRARY* constants in the file with values according to your current CUDA SDK installation.
### Installing
Once *install.py* has built all the needed files, you can copy the whole package to a place included in your PYTHONPATH (or change it to add this package's parent folder). See [this](http://superuser.com/questions/247620/how-to-globally-modify-the-default-pythonpath-sys-path) if you have any problem modifying it.
### Testing
Once installed you can run the tests in *pyRMSD/test* using:  
  
    \> python -m unittest testCondensedMatrix testMatrixHandler testMatrixNeighbours testMatrixStatistics testRMSDCalculators testPdbReader
  
You can use the benchmarks in *pyRMSD/benchmark* to tune CUDA or OpenMP parameters.  
## USING IT
### Getting coordinates
To use the module the first thing will be to extract all the coordinates from a PDB file. Coordinates must be stored in numpy arrays, using the same layout that  Prody uses:  

    Coordset: [Conformation 1, Conformation 2, ..., Conformation N]  
    Conformation: [Atom 1, Atom 2,..., Atom M]  
    Atom: [x,y,z]  

In order to do this there's a convenience class function in *pyRMSD/utils/proteinReading.py* called **Reader**. This will read a pdb file using the built in reader ("LITE_READER") or Prody ("PRODY_READER"), allowing to use its powerful selection language.  
  
    from pyRMSD.utils.proteinReading import Reader    
    reader = Reader("LITE_READER").readThisFile("my_trajectory.pdb").gettingOnlyCAs()
    coordinates = lite_reader.read()
    num_of_atoms = reader.numberOfAtoms
    num_of_framesreader.numberOfFrames
    
  
See 'pyRMSD/pyRMSD/test/testPdbReader.py for a simple usage example.
### Calculating the RMSD matrix
To calculate the RMSD matrix you can use directly the RMSD calculation function or a **MatrixHandler**or using directly one calculator object to feed a CondensedMatrix. A RMSDCalculator not only allows to create matrices, but also:  
* Pairwise RMSD calculation
* Reference vs the rest of the set.
* Reference vs following conformations.

This is a code snippet that shows one possible use scenario in which the calculator generates the all vs all rmsd matrix (with ALL the pairwise superpositions):
  
    import pyRMSD.RMSDCalculator
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets_mini, "PYTHON_CALCULATOR")
    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix = CondensedMatrix(rmsd)

As the resulting matrix is symmetric and its diagonal is 0, the rmsd_matrix object will store only the upper diagonal triangle (condensed matrix), in the same way [scipy.spatial.distance.pdist](http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
does.  
The calculator types can be one of these:  

* KABSCH_PYTHON_CALCULATOR
* QTRFIT_SERIAL_CALCULATOR
* QTRFIT_OMP_CALCULATOR
* QCP_CUDA_CALCULATOR (in CUDA capable machines)
* QCP_SERIAL_CALCULATOR
* QCP_OMP_CALCULATOR

Programatically, available calculators can be seen with:  
    
    from pyRMSD.availableCalculators import availableCalculators
    print availableCalculators()

### Matrix handlers
A **MatrixHandler** object will help you to create the matrix and will also help you saving and loading matrix data to disk.  

    from pyRMSD.matrixHandler import MatrixHandler  
    mHandler = MatrixHandler()  
    # Create a matrix with the coordsets and using a calculator  
    matrix = mHandler.createMatrixWithReader("amber_35k.pdb", "LITE_READER", "QCP_CUDA_CALCULATOR")  
    # Save the matrix to 'to_this_file.bin'  
    m_handler.saveMatrix("to_this_file")  
    # Load it from 'from_this_file.bin'  
    m_handler.loadMatrix("from_this_file")  
    # Get the inner CondensedMatrix instance
    rmsd_matrix = m_handler.getMatrix()  

### Accessing the RMSD matrix
You can access a matrix object contents like this:  

    rmsd_at_pos_2_3 = rmsd_matrix[2,3]

The row_lenght parameter will give you the... row length. Remember that the matrix is square and symmetric, so row_length == column_length, rmsd_matrix[i,j] == rmsd_matrix[j,i] and as it is a distance matrix, rmsd_matrix[i,i] == 0.

### Matrix statistics
The CondensedMatrix class also offers an efficient way to ask for the most common statistical moments. Use the methods **calculateMean**, **calculateVariance**, **calculateSkewness** and **calculateKurtosis** to get mean, variance, skewness and kurtosis ( easy, isn't it :) ). You can also use **calculateMax** and **calculateMin** to get the maximum and minimum value of the matrix.

## FUTURE IMPROVEMENTS
If you have used this package and you feel something is missing/incorrect or whatever, you can change it and contribute. Some examples of things that need to be improved are:  
* Solving bug in the CondensedMatrix object (erroneous creation when using a numpy array)
* Adding number of threads option for any OpenMP calculator.  **DONE**
* Adding  number of blocks and threads per block option in CUDA calculator.  **DONE**
* Create an installer using Python distutils (difficult because of the use of CUDA).  
* Add more tests.  
* Add more comments...  
* and improving this README!!  

##CREDITS
- Some Numpy helper functions were first seen in  http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays, by Lou Pecora (if I'm not wrong).

- The Python implementation of superposition was extracted from Prody source code (by [Ahmet Bakan](http://www.csb.pitt.edu/People/abakan/)) and modified, with the only goal of providing a python example to compare performance and stability.

- QCP superposition method code was adapted from the code [here](http://theobald.brandeis.edu/qcp/)

- The statistics function code was adapted from the work of jjhaag@dreamincode.net (available [here](http://www.dreamincode.net/code/snippet1447.htm) ).
