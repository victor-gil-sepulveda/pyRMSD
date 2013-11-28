pyRMSD
======

pyRMSD goal is the fast (and easy!) calculation of rmsd collective
operations, specially matrices of large ensembles of protein
conformations. It also offers a symmetric distance matrix implementation
with improved access speed and memory efficiency.

If you like it and you are using it for your scientific projects, please
cite `the pyRMSD
paper <http://bioinformatics.oxfordjournals.org/content/29/18/2363>`_.

pyRMSD distributed under MIT license, and it is currently on its version
4.0 .

Summary
=======

-  `1 - Features <#1---features>`_

   -  `Collective operations <#collective-operations>`_
   -  `Condensed matrix <#condensed-matrix>`_

-  `2 - Usage <#2---usage>`_

   -  `Getting coordinates <#getting-coordinates>`_
   -  `Calculating the RMSD matrix <#calculating-the-rmsd-matrix>`_
   -  `Available calculators <#available-calculators>`_
   -  `Matrix handlers <#matrix-handlers>`_
   -  `Accessing the RMSD matrix <#accessing-the-rmsd-matrix>`_
   -  `Matrix statistics <#matrix-statistics>`_

-  `3 - Building & Installation <#3---building--installation>`_

   -  `Before installation <#before-installation>`_
   -  `Linux and MacOs <#linux-and-macos>`_
   -  `Windows <#windows>`_

-  `4 - The custom building script <#4---the-custom-building-script>`_

   -  `Unix-based systems <#unix-based-systems>`_
   -  `Windows systems <#windows-systems>`_

      -  `Modifying system variables <#modifying-system-variables>`_

-  `5 - Testing (Developers) <#5---testing-developers>`_
-  `6 - Benchmarks (Developers) <#6---benchmarks-developers>`_
-  `Future improvements <#future-improvements>`_
-  `Credits <#credits>`_

1 - Features
------------

Collective operations
~~~~~~~~~~~~~~~~~~~~~

pyRMSD currently has 5 basic operations:
 1 - Pairwise RMSD calculation
 2 - One vs. following (of a sequence of conformers).
 3 - One vs. all the other conformations (of a sequence of conformers).
 4 - Pairwise RMSD matrix
 5 - Iterative superposition of a sequence.

All methods can use the same coordinates for fitting and RMSD
calculation, or a different set of coordinates for fitting (superposing)
and calculating RMSD (referred into the code as 'calculation
coordinates' ).

Currently pyRMSD implements a total of 3 superposition algorithms
(Kabsch's,QTRFIT and QCP) which can have serial or parallel versions
(OpenMP and CUDA in one case).

The available calculators so far are: \* KABSCH\_SERIAL\_CALCULATOR \*
KABSCH\_OMP\_CALCULATOR \* QTRFIT\_SERIAL\_CALCULATOR \*
QTRFIT\_OMP\_CALCULATOR \* QCP\_SERIAL\_CALCULATOR \*
QCP\_OMP\_CALCULATOR \* QCP\_CUDA\_CALCULATOR (in CUDA capable
machines*)* QCP\_CUDA\_MEM\_CALCULATOR (in CUDA capable machines\*)

In addition it offers 2 other calculators that do not perform
superposition (for cases in which the parts of interest of the system
are already superposed): \* NOSUP\_SERIAL\_CALCULATOR
\* NOSUP\_OMP\_CALCULATOR
This calculator will also center the coordinates, adding a little
unnecessary overhead. This overhead will be totally diluted when
calculating RMSD matrices though.

Finally it also holds a hidden calculator,
QCP\_SERIAL\_FLOAT\_CALCULATOR, maninly used to test against
QCP\_CUDA\_CALCULATOR in its float version.

Methods 1, 2 and 3 can be used to modify the input coordinates (the
input coordinates will be superposed). The iterative superposition
method will always have this behaviour as it would be senseless
otherwise. Conversely, RMSD matrix will never modify input coordinates.

pyRMSD can also have fitting symmetries and rotational calculation
symmetries into account. Documentation about this is on its way.

If you think you need new features to be added (or better examples)
click `here <#contact_features>`_.
\* Computing capability of the GPU must be equal or higher than 1.1
(>1.2 if built with double precision support).
###Condensed matrix
pyRMSD contains also a C written data type called CondensedMatrix. This
is a representation of a squared symmetric matrix and it will save you
half of the, otherwise redundant, memory. Besides, its write and read
access outperforms other implementations like pure python's list-based
and even Cython implementations (see the benchmarks folder). This means
that it will speed up for free any application that heavily relies on
accessing a distance matrix, like clustering algorithms. See the
examples below to get more insight about how to use it. ##2 - Usage
Some code snippets and explanations about them will be shown below. Note
that, as the code changes rapidly, this snippets can be outdated. I will
put all my effort for this not to happen, but if you detect that the
code examples are being more problematic than helpful for you, please
`contact me <#contact_features>`_. You will also find method and
variables documentation in the code. Do not hesitate to ask for more
documentation if you find is missing.

Getting coordinates
~~~~~~~~~~~~~~~~~~~

To use the module the first thing will be to extract all the coordinates
from a PDB file. Coordinates must be stored in numpy arrays, using the
same layout that Prody uses:

::

    Coordset: [Conformation 1, Conformation 2, ..., Conformation N]  
    Conformation: [Atom 1, Atom 2,..., Atom M]  
    Atom: [x,y,z]  

In order to do this there's a convenience class function in
*pyRMSD/utils/proteinReading.py* called **Reader**. This will read a pdb
file using the built in reader. Despite this, we encourage the use of
`Prody <http://www.csb.pitt.edu/prody/>`_ if you need to do any kind of
selection/manipulation.

::

    from pyRMSD.utils.proteinReading import Reader    
    reader = Reader().readThisFile("my_trajectory.pdb").gettingOnlyCAs()
    coordinates = reader.read()
    num_of_atoms = reader.numberOfAtoms
    num_of_frames = reader.numberOfFrames

See 'pyRMSD/pyRMSD/test/testPdbReader.py for a simple usage example.
###Calculating the RMSD matrix
To calculate the RMSD matrix you can use a **MatrixHandler** or use
directly one calculator object to feed a CondensedMatrix.

Using **MatrixHandler** to get the RMSD pairwise matrix (given that we
already have read the coordinates) will look like this:

::

    from pyRMSD.matrixHandler import MatrixHandler
    rmsd_matrix = MatrixHandler()\
       .createMatrix(coordinates, 'QCP_OMP_CALCULATOR')

Calculating the matrix using directly the RMSDCalculator is a little bit
more verbose:

::

    import pyRMSD.RMSDCalculator
    calculator = pyRMSD.RMSDCalculator.\
                    RMSDCalculator(coordsets,\
                    "QCP_SERIAL_CALCULATOR")
    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix = CondensedMatrix(rmsd)

As the resulting matrix is symmetric and its diagonal is 0, the
rmsd\_matrix object will store only the upper diagonal triangle
(condensed matrix), in the same way
`scipy.spatial.distance.pdist <http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html>`_
does.
###Available calculators Programatically, available calculators can be
queried with:

::

    from pyRMSD.availableCalculators import availableCalculators
    print availableCalculators()

Matrix handlers
~~~~~~~~~~~~~~~

A **MatrixHandler** object will help you to create the matrix and will
also help you saving and loading matrix data to disk.

::

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

Accessing the RMSD matrix
~~~~~~~~~~~~~~~~~~~~~~~~~

You can access a matrix object contents like this:

::

    rmsd_at_pos_2_3 = rmsd_matrix[2,3]

The **row\_lenght** parameter will give you the... row length. Remember
that the matrix is square and symmetric, so row\_length ==
column\_length, rmsd\_matrix[i,j] == rmsd\_matrix[j,i] and as it is a
distance matrix, rmsd\_matrix[i,i] == 0.

Matrix statistics
~~~~~~~~~~~~~~~~~

The CondensedMatrix class also offers an efficient way to ask for the
most common statistical moments. Use the methods **calculateMean**,
**calculateVariance**, **calculateSkewness** and **calculateKurtosis**
to get mean, variance, skewness and kurtosis ( easy, isn't it :) ). You
can also use **calculateMax** and **calculateMin** to get the maximum
and minimum value of the matrix.

3 - Building & Installation
---------------------------

Before installation
~~~~~~~~~~~~~~~~~~~

**Users** only need to install Python version 2.6/2.7 (pyRMSD has only
been tested with those, however it may work with another versions of the
Python 2.X family). Numpy is also required. Surely you already have it
into your machine, but, in the case you don't, it can be found
`here <http://sourceforge.net/projects/numpy/files/>`_. There you will
be able to find installers for almost all the combinations of platforms
and Python versions you can think about.

**Developers** may remember that header files of Python and Numpy may be
accessible, and your Python installation must contain the python shared
library. This usually means that you have to compile it using
./configure --enable-shared before building Python (usually 2.7
distributions already come with this library). Prody is not a
dependency, but I encourage its use to handle coordinates, as it is
well-tested and powerful tool.

Linux and MacOs
~~~~~~~~~~~~~~~

Those users have the following choices:

**1)** Using the 'setup.py' file inside the root folder by typing:

::

    > python setup.py install

with superuser provileges. Use 'build' instead to only build it. This is
the usual way python packages are deployed. As 'distutils' do not
support CUDA directly, your package will not be able to use CUDA
calculators.

**2)** Using the custom build.py script in pyRMSD main folder. This will
compile a version of pyRMSD following your configuration details. To
finish the installation you will need to change your PYTHONPATH in order
to point to the parent folder of the package (or copy it in a folder
already inside your PYTHONPATH). See
`this <http://superuser.com/questions/247620/how-to-globally-modify-the-default-pythonpath-sys-path>`_
if you have any problem modifying it.

Windows
~~~~~~~

Windows Installation is discontinued. I keep some very basic
instructions `here <#windows-systems>`_ though.

4 - The custom building script
------------------------------

pyRMSD includes a small build script that is indeed a recipe to compile
the C extensions of pyRMSD. The build.py script is the most versatile
way to compile pyRMSD and will work in almost all situations. With this
script one can build x86 and x64 distributions with enabled or disbled
CUDA calculators. Invoke it from pyRMSD root folder with:

::

    > python build.py \[OPTIONS\]

By default this script won't do anything. OPTIONS can be one of these:

--build -> to compile pyRMSD (OpenMP version).
--cuda single/double -> to compile it with single or double precission
(you must specify only one). Double precission will not work in old
cards even if they are CUDA capable.
--clean -> Will remove any generated .o files. --build --clean is a good
combination if you are not a developer.
--clean-all -> Will remove all generated files. Combine this one with
any other is not a good idea. It will remove any useful built file.
--build-conf -> Will determine the file (inside *build\_conf* folder)
storing the configuration info.
--help/-h -> Will write some hints about the options.

This script uses distutil's *sysconfig* package to get the search path
for python headers and libs automatically.
The building script will try to guess the location of the needed files
for compilation, but it can be easily modified to be able to handle all
kind of scenarios.

Configuration files
~~~~~~~~~~~~~~~~~~~

As stated before, multiple configuration files can be used by the
building script to feed it with the correct variables. This
configuration files are stored in the *build\_conf* folder and by
default, the 'default.conf' file is loaded (equivalent to *--build-conf
default.conf*). These are the parameters that can be changed. If one key
is not present in the file, then the contents of the key inside the
'default.conf' file are used.

-  "CUDA\_BASE": Base folder of the CUDA distribution installed.
-  "CUDA\_INCLUDE\_FOLDER": Folder inside CUDA\_BASE where CUDA headers
   are found.
-  "CUDA\_LIBRARIES\_FOLDER": Folder inside CUDA\_BASE where CUDA
   libraries are found.
-  "CUDA\_ARCHITECHTURE": Arquitecture of the GPU (in the sm\_XX
   format).
-  "CUDA\_LIBRARY": Name of the cuda library (usually 'cudart').
-  "CUDA\_OPTIONS": Some options for the CUDA compiler.
-  "DEFINE\_USE\_CUDA": Allows to redefine the macro that tells the
   preprocessor that CUDA is going to be used. Changing this means you
   also changed parts of the code, so is not adviced.
-  "PYTHON\_EXTENSION\_OPTIONS": Compiler options usually added to
   compile Python extensions.
-  "PYTHON\_INCLUDE\_FOLDER": If "AUTO" it will use
   'distutils.sysconfig' to obtain the location of Python's header
   files, if not it must be the full location of python's header files
   folder.
-  "PYTHON\_LIBRARY\_FOLDER": If "AUTO" it will use
   'distutils.sysconfig' to obtain the location of Python's libraries.
   If "AUTO\_ALT", it will use os.path.dirname (useful for those
   situations in which 'distutils.sysconfig' fails to get the propper
   folder). The other option is the full location of python's libraries
   folder.
-  "PYTHON\_LIBRARY" : The python library to link against.
-  "OPENMP\_OPTION" : The openmp flag for the compiler used (nowadays
   the only compiler is gcc, so this flag must be -fopenmp in almost all
   cases).
-  "NUMPY\_INCLUDE" : If "AUTO" it will use ''numpy.get\_include' to get
   the folder where numpy headers are. If not it must be the full path
   of this folder.
-  "PYTHON\_EXTENSION\_LINKING\_OPTIONS": Compiler options usually added
   to link Python extensions.

See the examples into the *build\_conf* if you want to create your own
configuration files.

Unix-based systems
~~~~~~~~~~~~~~~~~~

The building script was used in a Ubuntu x86 and Ubuntu x64 Os, as well
as in MacOs (Snow Leopard) to perform a non CUDA build. It had the
*python-dev* package installed, so python headers were available.
PYTHON\_X constants were left unchanged.
It was also used under Ubuntu x64 with CUDA 4.2 to build the CUDA
enabled version.
####Mac Users Roman Sloutsky warns that if you're not able to compile
using the build script with default configuration options, just try to
change "PYTHON\_LIBRARY\_FOLDER":"AUTO" to
"PYTHON\_LIBRARY\_FOLDER":"AUTO\_ALT" in "default.conf". Creating a new
configuration file with only this entry will also work.

Windows systems
~~~~~~~~~~~~~~~

A preliminary version of the build script was also tested in Windows 7
32 and 64 systems using MinGW compiler tools. Here are the steps
followed to succesfully compile the extensions:

- `Download <http://www.mingw.org/>`_ and install MinGW. Then add its
/bin folder to Windows PATH
- `Download <http://www.python.org/download/releases/2.7.3/>`_ and
install Python 2.7.3
- `Download <http://www.scipy.org/Download>`_ and install Numpy (tested
with v. 1.7.0 for python 2.7)
- [Optional] `Download <http://www.csb.pitt.edu/prody/getprody.html>`_
and install Prody (tested with v. 1.4.1 for python 2.7)

*PYTHON\_INCLUDE\_FOLDER* and *PYTHON\_LIBRARY\_FOLDER* constants were
changed to match our Python installation paths.
*PYTHON\_EXTENSION\_LINKING\_OPTIONS* and *PYTHON\_EXTENSION\_OPTIONS*
were also changed to fit Windows extension creation options.
Once everything is built, create (or modify) the PYTHONPATH system
variable and make it point the pyRMSD folder.

Modifying system variables
^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to create or modify a system variable under Windows 7, you will
have to go to Control Panel -> System and Security -> System -> Advanced
System Settings.

5 - Testing (Developers)
------------------------

Once installed you can run the tests in *pyRMSD/test* using:

::

    > python -m unittest discover

Currently only the *test\_create\_with\_reader* test will fail if all
the dependencies are fullfilled (it's unwritten yet). If you didn't
build pyRMSD with CUDA support, 5 tests will be skipped.

If you compiled the package using the build script, an extra test suite
will be available in the *src/calculators/test* folder with pure C
tests. Run it inside this folder with ./test\_rmsdtools\_main

6 - Benchmarks (Developers)
---------------------------

Also available to users, inside the */benchmark* folder, there are the
benchmarks used to assest the performance of pyRMSD.
There one can find some small scripts to test OpenMP parametrizations,
calculation time of every implementation or even a small floating point
error check.

Future improvements
-------------------

If you have used this package and you feel something is
missing/incorrect or whatever, you can change it and contribute. Some
examples of things that need to be improved are:
\* Solving bug in the CondensedMatrix object (erroneous creation when
using a numpy array)
\* Adding number of threads option for any OpenMP calculator. **DONE**
\* Adding number of blocks and threads per block option in CUDA
calculator. **DONE**
\* Create a unique installer using Python distutils (difficult because
of the use of CUDA).
\* C code needs more comments and to encapsulate function arguments (an
absolutely needed major refactoring). **DONE** \* Matrix generation load
balance can be improved. \* Names in C code must be more
self-explanatory. **IMPROVED** \* C code must be revised and further
simplified. **DONE**
\* Add more and better tests to increase coverage.
\* Convert build.py in a Makefile generator. **PARTIALLY DONE** (it now
acts just like make would do)
\* Add more comments...
\* Symmetry features need to be explained. \* and of course improving
this README!! If you want to add new features (for instance mass
weighting) do not hesitate to contact me at: victor.gil.sepulveda at
gmail.com. Of course, you can fork the repository and add as many
features and improvements as you want.

Credits
-------

-  Some Numpy helper functions were first seen in
   http://www.scipy.org/Cookbook/C\_Extensions/NumPy\_arrays, by Lou
   Pecora (if I'm not wrong).

-  The initial Python implementation of superposition was extracted from
   Prody source code (by `Ahmet
   Bakan <http://www.csb.pitt.edu/People/abakan/>`_) and modified, with
   the only goal of providing a python example to compare performance
   and stability. The iterative superposition algorithm is a direct
   translation of his iterpose algorithm.

-  QCP superposition method code was adapted from the code
   `here <http://theobald.brandeis.edu/qcp/>`_

-  The statistics function code was adapted from the work of
   jjhaag@dreamincode.net (available
   `here <http://www.dreamincode.net/code/snippet1447.htm>`_ ).

-  Kabsch algorithm code was adapted from the work of `Dr. Bosco K.
   Ho <http://boscoh.com/>`_. I would like to give him special tanks for
   his help.

-  As far as I know the first CUDA implementation of QCP is from `Yutong
   Zhao <http://proteneer.com/blog/>`_. He went a step further trying to
   improve memory coalescence by changing the coordinates overlay.
   Pitifully his code is not open source.


