pyRMSD
==================
pyRMSD goal is the fast (and easy!) calculation of rmsd matrices of large ensembles of protein conformations. It also offers a symmetric matrix 
representation which fosters memory access speed over access security and also is very memory efficient. 
The use of this matrix object can speed up any algorithm using a rmsd matrix a 4-5X. 

Building and Installing
==================
Dependencies
-------------
Before installing be sure you have Python 2.7, [numpy](http://numpy.scipy.org/) and [scipy](http://www.scipy.org/) (if you want to use the tests) and [prody](http://pypi.python.org/pypi/ProDy/) 
(if you want to be able parse PDBs... which will be surely the case). All this packages are really easy to install (well... scipy can be a little bit tricky in some
systems).
Building
-----------
To build the code (over Linux, no Windows support yet) you must execute *install.py*.  
    python install.py  
This will build all serial and OpenMP calculators. *install.py* can be modified to add/remove compiling options (and you know what you're doing).  
If you want to add the CUDA calculators, just use:  
    python install.py --cuda  
and it will also try to build the available CUDA calculators. In this case you will surely have to modify the *CUDA_BASE*, *CUDA_INCLUDE*, *CUDA_LIB*, *CUDA_ARCH* and *CUDA_LIBRARY* constants in the file with values according to your CUDA SDK installation.
Installing
------------
Once *install.py* has built all the needed files, you can copy the whole package to a place included in your PYTHONPATH (or change it to add this package's parent folder). See [this](http://superuser.com/questions/247620/how-to-globally-modify-the-default-pythonpath-sys-path) if you have problems modifying it.
Testing
--------
Once installed you can run the tests in *pyRMSD/test* using:  
    python -m unittest  

USING IT
=========
To use the module the first thing will be to extract all the coordinates from a PDB file. Coordinates must be stored in the same 
storage format that prody does. So if you don't have it yet, get it! (unless you want to replicate their job :P).
In fact, you can use prody through a wrapper in 'utils.py' (see 'pyRMSD/pyRMSD/test/test.py for a usage example).

FUTURE IMPROVEMENTS
============
If you have used this package and you feel something is missing/incorrect or whatever, you can change it and contribute. Some examples of things that need to be improved are:  
* Adding number of threads option for any OpenMP calculator.  
* Adding  number of blocks and threads per block option in CUDA calculator.  
* Create an installer using Python distutils (difficult because of the use of CUDA).  
* Add more tests.  
* Add more comments!!  

(coordinates_set, number_of_conformations, atoms_per_conformation)  = getCoorsetsFromPDB(pdb_path)
--------------------------------------------------------------------------------------------------
Will return a tuple with all the coordinates of all the frames inside the PDB, the number of conformations, and the number of 
atoms per conformation. 'pdb_path' is the path of the PDB file to be parsed.

This coordinates can directly feed the two actual functions you will need for rmsd calculation. This 2 functions are located 
into the 'RMSD.py' module and are:


rmsd = oneVsTheOthers(target,coordsets,calcType)
-----------------------------------------
Will calculate the RMSD of conformation 'target' (an integer from 0 to the number of conformations) with the conformations from 
'target'+1 to number_of_conformations. It will return a vector with the corresponding rmsds.


rmsd = calculateRMSDCondensedMatrix(coordsets,calcType)
-------------------------------------------------------
Will calculate the all vs all rmsd matrix (with ALL the pairwise superpositions). As the resulting matrix is symmetric and its 
diagonal is 0, it will store only the upper diagonal triangle (condensed matrix), in the same way scipy.spatial.distance.pdist
does (http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html).


For both functions possible calculator types (calcType) are:
"PYTHON_CALCULATOR" 	- To use the Python implementation.
"SERIAL_CALCULATOR" 	- To use the C Serial implementation.
"OMP_CALCULATOR" 	- To use the C Serial implementation with OpenMP.
"CUDA_CALCULATOR" 	- To use the CUDA implementation.



TESTING
========
In order to test your build go to the folder pyRMSD/pyRMSD/test, uncompress amber_long.tar.bz and run test.py (python test.py). It will 
test all the implementations and do a benchmark on your computer.


Credits and Thanks
==================
- Yutong Zhao for its cuda rmsd code (unreleased yet). One can get the compiled objets from his blog (http://proteneer.com/blog/) .

- Manuel Rivero, Ryoji Takahasi and Israel Cabeza de Vaca for the initial code for the serial rmsd code (part of PELE++ project).

- Helper functions where extracted from http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays, by Lou Pecora if I'm not wrong.

- The Python implementation of superposition was extracted from Prody source code (by Ahmet Bakan [http://www.csb.pitt.edu/People/abakan/]), with
the only goal of providing a framework for improving the code.

- All the other stuff (wrapping, improvements, etc) was done by myself: Víctor Gil Sepúlveda.

I hope it will be helpful for you all.
