# pyRMSD
pyRMSD goal is the fast (and easy!) calculation of rmsd collective operations, specially matrices of large ensembles of protein conformations. It also offers a symmetric distance matrix implementation with improved access speed and memory efficiency.

pyRMSD distributed under MIT license, and it is currently on its version 3.0 .

- [1 - Features](#1---features)  
- [2 - Usage](#2---usage)  
	- [Getting coordinates](#getting-coordinates)  
	- [Calculating the RMSD matrix](#calculating-the-rmsd-matrix)  
	- [Available calculators](#available-calculators)  
	- [Matrix handlers](#matrix-handlers)  
	- [Accessing the RMSD matrix](#accessing-the-rmsd-matrix)
	- [Matrix statistics](#matrix-statistics)  
- [3 - Building & Installation](#3---building--installation)  
	- [Before installation](#before-installation)  
	- [Linux and MacOs](#linux-and-macos)  
	- [Windows](#windows)  
- [4 - The custom building script](#4---the-custom-building-script)  
	- [Unix-based systems](#unix-based-systems)  
	- [Windows systems](#windows-systems)  
		- [Modifying system variables](#modifying-system-variables)  
- [5 - Testing (Developers)](#5---testing-developers)  
- [6 - Benchmarks (Developers)](#6---benchmarks-developers)  
- [Future improvements](#future-improvements)  
- [Credits](#credits)  
    
##1 - Features  
pyRMSD currently has 5 basic operations:  
1 - Pairwise RMSD calculation  
2 - One vs. following (in a sequence of conformers).  
3 - One vs. all the other conformations (in a sequence of conformers).  
4 - Pairwise RMSD matrix  
5 - Iterative superposition of a sequence.  

All methods can use the same coordinates for fitting and RMSD calculation, or a different set of coordinates for fitting (superposing) and calculating RMSD.

In addition, methods 1, 2 and 3 can be used to modify the input coordinates (the input coordinates will be superposed). The iterative superposition method will always have this behaviour as it would be senseless otherwise. 
  
##2 - Usage  
Some code snippets and explanations about them will be shown below. Note that as the code changes rapidly, this snippets can be outdated. I will put all our effort for this not to happen, but if you detect that anything is not actually working, please contact me.
  
###Getting coordinates
To use the module the first thing will be to extract all the coordinates from a PDB file. Coordinates must be stored in numpy arrays, using the same layout that Prody uses:  

    Coordset: [Conformation 1, Conformation 2, ..., Conformation N]  
    Conformation: [Atom 1, Atom 2,..., Atom M]  
    Atom: [x,y,z]  

In order to do this there's a convenience class function in *pyRMSD/utils/proteinReading.py* called **Reader**. This will read a pdb file using the built in reader.  
  
    from pyRMSD.utils.proteinReading import Reader    
    reader = Reader().readThisFile("my_trajectory.pdb").gettingOnlyCAs()
    coordinates = reader.read()
    num_of_atoms = reader.numberOfAtoms
    num_of_framesreader.numberOfFrames

See 'pyRMSD/pyRMSD/test/testPdbReader.py for a simple usage example.
###Calculating the RMSD matrix
To calculate the RMSD matrix you can use a **MatrixHandler** or use directly one calculator object to feed a CondensedMatrix. 

Using **MatrixHandler** to get the RMSD pairwise matrix (given that we already have read the coordinates) will look like this:
  
    from pyRMSD.matrixHandler import MatrixHandler
    rmsd_matrix = MatrixHandler()\
       .createMatrix(coordinates, 'QCP_OMP_CALCULATOR')
  
Calculating the matrix using directly the RMSDCalculator is a little bit more verbose:

    import pyRMSD.RMSDCalculator
    calculator = pyRMSD.RMSDCalculator.\
                    RMSDCalculator(coordsets,\
                    "QCP_SERIAL_CALCULATOR")
    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix = CondensedMatrix(rmsd)

As the resulting matrix is symmetric and its diagonal is 0, the rmsd_matrix object will store only the upper diagonal triangle (condensed matrix), in the same way [scipy.spatial.distance.pdist](http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
does.  
###Available calculators
The calculator type can be one of these:  

* KABSCH_SERIAL_CALCULATOR
* KABSCH_OMP_CALCULATOR
* QTRFIT_SERIAL_CALCULATOR
* QTRFIT_OMP_CALCULATOR
* QCP_SERIAL_CALCULATOR
* QCP_OMP_CALCULATOR
* QCP_CUDA_CALCULATOR (in CUDA capable machines *)
* QCP_CUDA_MEM_CALCULATOR (in CUDA capable machines *)

Which implement Kabsch's superposition algorithm, QTRFIT, and Theobald's QCP.

Programatically, available calculators can be seen with:  
    
    from pyRMSD.availableCalculators import availableCalculators
    print availableCalculators()
  
\* Computing capability of the GPU must be equal or higher than 1.1 (>1.2 if built with double precision support).

###Matrix handlers
A **MatrixHandler** object will help you to create the matrix and will also help you saving and loading matrix data to disk.   

    from pyRMSD.matrixHandler import MatrixHandler  
    # Create a matrix with the coordsets and using a calculator
    mHandler = MatrixHandler()  
    matrix = mHandler.createMatrix( coordsets,\
                                    "QCP_CUDA_CALCULATOR")
                                    
    # Save the matrix to 'to_this_file.bin'  
    m_handler.saveMatrix("to_this_file")  
    
    # Load it from 'from_this_file.bin'  
    m_handler.loadMatrix("from_this_file") 
    
    # Get the inner CondensedMatrix instance
    rmsd_matrix = m_handler.getMatrix()  

###Accessing the RMSD matrix
You can access a matrix object contents like this:  

    rmsd_at_pos_2_3 = rmsd_matrix[2,3]

The **row_lenght** parameter will give you the... row length. Remember that the matrix is square and symmetric, so row_length == column_length, rmsd_matrix[i,j] == rmsd_matrix[j,i] and as it is a distance matrix, rmsd_matrix[i,i] == 0.

###Matrix statistics
The CondensedMatrix class also offers an efficient way to ask for the most common statistical moments. Use the methods **calculateMean**, **calculateVariance**, **calculateSkewness** and **calculateKurtosis** to get mean, variance, skewness and kurtosis ( easy, isn't it :) ). You can also use **calculateMax** and **calculateMin** to get the maximum and minimum value of the matrix.

##3 - Building & Installation
###Before installation
**Users** only need to install Python version 2.6/2.7 (pyRMSD has only been tested with those, however it may work with another versions of the Python 2.X family). Numpy is also required. Surely you already have it installed in your machine, but in the case you didn't it can be found [here](http://sourceforge.net/projects/numpy/files/), where you will be able to find installers for almost all the combinations of platforms and Python versions you can think about.

**Developers** may remember that header files of Python and Numpy may be accessible, and your Python installation must contain the python shared library. This usually means that you have to compile it using ./configure --enable-shared before building Python (usually 2.7 distributions already come with this library). Prody is not a dependency, but I encourage its use to handle coordinates, as it is well-tested and powerful tool.

###Linux and MacOs
  
Those users have the following choices:  
  
**1)** Using the 'setup.py' file inside the root folder by typing:  

    > python setup.py install
  
(or 'build') to only build it), which is the usual way python packages are deployed. As 'distutils' do not support CUDA directly, your package will not be able to use CUDA calculators.

**2)** Using the custom build.py script in pyRMSD main folder with:  
  
    > python build.py  
  
or  
  
    > python build.py --cuda single/double  
  
The build.py script is the most versatile way to compile pyRMSD and will work in almost all situations. With this script one can build x86 and x64 distributions with enabled CUDA calculators if the --cuda flag is used (followed by **single** or **double** depending on the precision of the floating point operations you want to use / your GPU allows). There is more information about the build.py script [here](#build_linux). 

###Windows  
  
Windows users have the following choices:  
  
**1)** Using a precompiled windows installer (in 'wininst/pyRMSD-3.0.win32-py2.7.msi'). Those are used as any other regular Windows Installer (double click the executable and follow instructions).  
  
**2)** Using the custom build.py script in pyRMSD main folder with:  
  
    > python build_windows.py
  
Please look [here](#build_win) if you need further iformation about the windows version of the custom build script.  

##4 - The custom building script  
pyRMSD includes a small build script that is indeed a recipe to compile the C extensions of pyRMSD. This script uses distutil's *sysconfig* package to get the search path for python headers and libs automatically.  
The building script will try to guess the location of the needed files for compilation, but it can be easily modified to be able to handle all kind of scenarios.  

###Unix-based systems  
The script was used in a Ubuntu x86 and Ubuntu x64 Os, as well as in MacOs (Snow Leopard) to perform a non CUDA build. It had the *python-dev* package installed, so python headers were available. PYTHON_X constants were left unchanged.  
It was also used under Ubuntu x64 with CUDA 4.2 to build the CUDA enabled version.  
If you are going to use it to build a CUDA enabled version you may have to change the *CUDA_BASE* constant, which needs to point to the base directory of your CUDA installation (in our case  */usr/local/cuda-4.2*). Required headers and libs are usually stored inside the */include* and */lib64* (*/lib* in x86 systems) subfolders, but you can also change it by modifying *CUDA_INCLUDE_FOLDER* and *CUDA_LIBRARIES_FOLDER*. Change *CUDA_ARCHITECHTURE* to match the architecture of your GPU.  
Finally you will need to change your PYTHONPATH in order to point to the parent folder of the package (or copy it in a folder already inside your PYTHONPATH). See [this](http://superuser.com/questions/247620/how-to-globally-modify-the-default-pythonpath-sys-path) if you have any problem modifying it.  

###<a id="build_win"></a>Windows systems  
The build script has also been tested in Windows 7 32 and 64 systems using MinGW compiling tools. Here are the steps followed to succesfully compile the extensions:  
  
\- [Download](http://www.mingw.org/) and install MinGW. Then add its /bin folder to Windows PATH  
\- [Download](http://www.python.org/download/releases/2.7.3/) and install Python 2.7.3  
\- [Download](http://www.scipy.org/Download) and install Numpy (tested with v. 1.7.0 for python 2.7)    
\- \[Optional\] [Download](http://www.csb.pitt.edu/prody/getprody.html) and install Prody (tested with v. 1.4.1 for python 2.7)
  
Inside the build_windows.py, *PYTHON_INCLUDE_FOLDER* and *PYTHON_LIBRARY_FOLDER* constants were changed to match our Python installation paths.  
*PYTHON_EXTENSION_LINKING_OPTIONS* and *PYTHON_EXTENSION_OPTIONS* were also changed to fit Windows extension creation options.    
Once everything is built, create (or modify) the PYTHONPATH system variable and make it point the pyRMSD folder.  

#### Modifying system variables  

In order to create or modify a system variable under Windows 7, you will have to go to Control Panel -> System and Security -> System -> Advanced System Settings.

##5 - Testing (Developers)  
Once installed you can run the tests in *pyRMSD/test* using:  
  
    > python -m unittest discover
  
Currently only the *test_create_with_reader* test will fail if all the dependencies are fullfilled (it's unwritten yet). 
If you didn't build pyRMSD with CUDA support, 5 tests will be skipped.

If you compiled the package using the build script, an extra test suite will be available in the src/calculators/test folder with pure C tests.  

##6 - Benchmarks (Developers)
Also available to users, inside the */benchmark* folder, there are the benchmarks used to assest the performance of pyRMSD.  
There one can find some small scripts to test OpenMP parametrizations, calculation time of every implementation or even a small floating point error check.
  
##Future improvements
If you have used this package and you feel something is missing/incorrect or whatever, you can change it and contribute. Some examples of things that need to be improved are:  
* Solving bug in the CondensedMatrix object (erroneous creation when using a numpy array)  
* Adding number of threads option for any OpenMP calculator.  **DONE**  
* Adding  number of blocks and threads per block option in CUDA calculator.  **DONE**  
* Create an installer using Python distutils (difficult because of the use of CUDA).  
* C code needs more comments and to encapsulate funtion arguments.  
* Add more tests.  
* Add more comments...  
* and improving this README!!  

##Credits
- Some Numpy helper functions were first seen in  http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays, by Lou Pecora (if I'm not wrong).

- The initial Python implementation of superposition was extracted from Prody source code (by [Ahmet Bakan](http://www.csb.pitt.edu/People/abakan/)) and modified, with the only goal of providing a python example to compare performance and stability. The iterative superposition algorithm is a direct translation of his iterpose algorithm.

- QCP superposition method code was adapted from the code [here](http://theobald.brandeis.edu/qcp/)

- The statistics function code was adapted from the work of jjhaag@dreamincode.net (available [here](http://www.dreamincode.net/code/snippet1447.htm) ).

- Kabsch algorithm code was adapted from the work of [Dr. Bosco K. Ho](http://boscoh.com/)  
- 
