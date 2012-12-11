# pyRMSD
pyRMSD goal is the fast (and easy!) calculation of rmsd collective operations, specially matrices of large ensembles of protein conformations. It also offers a symmetric matrix representation which fosters access speed over access security and also is very memory efficient.

## Building and Installing
### Dependencies

Before installing be sure you have Python 2.7 and [numpy](http://numpy.scipy.org/). Some soft dependencies are [scipy](http://www.scipy.org/), only necessary to execute one test (statistics), and [prody](http://pypi.python.org/pypi/ProDy/) which will add more power to pdb parsing (and also mandatory for one or two tests). All this packages are really easy to install (well... installing scipy can be a little bit tricky in some systems).

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
  
    > python -m unittest testCondensedMatrix testMatrixHandler testMatrixNeighbours testMatrixStatistics testRMSDCalculators testPdbReader
  
You can use the benchmarks in *pyRMSD/benchmark* to tune CUDA or OpenMP parameters.  
## USING IT
### Getting coordinates
To use the module the first thing will be to extract all the coordinates from a PDB file. Coordinates must be stored the same layout that  prody uses:  

    Coordset: [Conformation 1, Conformation 2, ..., Conformation N]  
    Conformation: [Atom 1, Atom 2,..., Atom M]  
    Atom: [x,y,z]  

In order to do this there's a convenience class function in *pyRMSD/utils/proteinReading.py* called **Reader**. This will read a pdb file using the built in reader ("LITE_READER") or Prody ("PPRODY_READER"), allowing to use its powerful selection language.  
  
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
The available calculator types in a CUDA capable machine can be one of these:  

* KABSCH_PYTHON_CALCULATOR
* QTRFIT_SERIAL_CALCULATOR
* QTRFIT_OMP_CALCULATOR
* QCP_CUDA_CALCULATOR
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
* Adding number of threads option for any OpenMP calculator.  **DONE**
* Adding  number of blocks and threads per block option in CUDA calculator.  **DONE**
* Create an installer using Python distutils (difficult because of the use of CUDA).  
* Add more tests.  
* Add more comments...  
* and improving this README!!  

##CREDITS
- Some Numpy helper functions were first seen in  http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays, by Lou Pecora (if I'm not wrong).

- The Python implementation of superposition was extracted from Prody source code (by [Ahmet Bakan](http://www.csb.pitt.edu/People/abakan/)) and modified, with the only goal of provide a python example to compare performance and stability.

- QCP superposition method code was adapted from the code [here](http://theobald.brandeis.edu/qcp/)

- The statistics function code was adapted from the work of jjhaag@dreamincode.net (available [here](http://www.dreamincode.net/code/snippet1447.htm).

