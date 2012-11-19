# pyRMSD
pyRMSD goal is the fast (and easy!) calculation of rmsd matrices of large ensembles of protein conformations. It also offers a symmetric matrix representation which fosters access speed over access security and also is very memory efficient.  
The use of this matrix object can speed up any algorithm using a rmsd matrix a 4-5X. 

## Building and Installing
### Dependencies

Before installing be sure you have Python 2.7, [numpy](http://numpy.scipy.org/) and [scipy](http://www.scipy.org/) (if you want to use the tests) and [prody](http://pypi.python.org/pypi/ProDy/) (if you want to be able parse PDBs... which will be surely the case). All this packages are really easy to install (well... scipy can be a little bit tricky in some
systems).
### Building
To build the code (over Linux, no Windows support yet) you must execute *install.py*.  
  
    > python install.py
  
This will build all serial and OpenMP calculators. *install.py* can be modified to add/remove compiling options (and you know what you're doing).  
If you want to add the CUDA calculators, just use:  
  
    >python install.py --cuda
  
and it will also try to build the available CUDA calculators. In this case you will surely have to modify the *CUDA_BASE*, *CUDA_INCLUDE*, *CUDA_LIB*, *CUDA_ARCH* and *CUDA_LIBRARY* constants in the file with values according to your CUDA SDK installation.
### Installing
Once *install.py* has built all the needed files, you can copy the whole package to a place included in your PYTHONPATH (or change it to add this package's parent folder). See [this](http://superuser.com/questions/247620/how-to-globally-modify-the-default-pythonpath-sys-path) if you have problems modifying it.
### Testing
Once installed you can run the tests in *pyRMSD/test* using:  
  
    > python -m unittest testCondensedMatrix testMatrixHandler testMatrixNeighbours testMatrixStatistics testRMSDCalculators  
  
## USING IT
### Getting coordinates
To use the module the first thing will be to extract all the coordinates from a PDB file. Coordinates must be stored the same layout that  prody uses:  
   Coordset: [Conformation 1, Conformation 2, ..., Conformation N]  
   Conformation: [Atom 1, Atom 2,..., Atom M]  
   Atom: [x,y,z]  
In order to do this there's a convenience function in *pyRMSD/utils/proteinReading.py* called **getCoordsetsFromPDB** (which will only pick **CA** atoms!!):  
    import pyRMSD.utils
    coordsets, number_of_conformations, number_of_atoms = pyRMSD.utils.proteinReading.getCoordsetsFromPDB('my.pdb',"a prody selection string")
See 'pyRMSD/pyRMSD/test/testRMSDCalculators.py for a simple usage example.
### Calculating the RMSD matrix
To calculate the RMSD matrix you can use directly the RMSD calculation function or a **MatrixHandler**.
Using the calculator is very straighforward. Use **calculateRMSDCondensedMatrix** in pyRMSD.RMSD to calculate the array representing a condensed matrix of the all vs all RMSDs, or **oneVsTheOthers** to calculate the rmsd values of conformation i vs conformations [i+1 to N] (This one is also very useful if you want to do a pairwise RMSD).  

    import pyRMSD.RMSD
    from pyRMSD.condensedMatrix import CondensedMatrix
    rmsd_array = pyRMSD.RMSD.calculateRMSDCondensedMatrix(coordsets, calculator_type)
    rmsd_matrix = CondensedMatrix(rmsd)

Will calculate the all vs all rmsd matrix (with ALL the pairwise superpositions). As the resulting matrix is symmetric and its 
diagonal is 0, it will store only the upper diagonal triangle (condensed matrix), in the same way scipy.spatial.distance.pdist
does (http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html). **calculator_type** can be one of theese:  

* PYTHON_CALCULATOR
* SERIAL_CALCULATOR
* OMP_CALCULATOR
* THEOBALD_CUDA_CALCULATOR
* THEOBALD_SERIAL_CALCULATOR
* THEOBALD_SERIAL_OMP_CALCULATOR

Available calculators can be seen with:
    
    import pyRMSD.RMSD
    print pyRMSD.RMSD.availableCalculators()

A **MatrixHandler** objecto 
### Accessing the RMSD matrix
You can access a matrix object contents like this:  
    rmsd_at_pos_2_3 = rmsd_matrix[2,3]
The row_lenght parameter will give you the... row length. Remember that the matrix is square and symmetric, so row_length == column_length, 
rmsd_matrix[i,j] == rmsd_matrix[j,i] and as it is a distance matrix, rmsd_matrix[i,i] == 0.

## FUTURE IMPROVEMENTS
If you have used this package and you feel something is missing/incorrect or whatever, you can change it and contribute. Some examples of things that need to be improved are:  
* Adding number of threads option for any OpenMP calculator.  
* Adding  number of blocks and threads per block option in CUDA calculator.  
* Create an installer using Python distutils (difficult because of the use of CUDA).  
* Add more tests.  
* Add more comments!!  

##CREDITS

- Helper functions where extracted from http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays, by Lou Pecora if I'm not wrong.

- The Python implementation of superposition was extracted from Prody source code (by [Ahmet Bakan](http://www.csb.pitt.edu/People/abakan/)), with the only goal of provide a python example to compare.

- QCP superposition method code was adapted from the code [here](http://theobald.brandeis.edu/qcp/)

- The statistics function code was adapted from the work of jjhaag@dreamincode.net (available [here](http://www.dreamincode.net/code/snippet1447.htm).

