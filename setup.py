'''
Created on 25/02/2013

@author: victor
'''
from distutils.core import setup, Extension
import numpy
setup(name='pyRMSD',
      version='2.0',
      description='pyRMSD is a small Python package that aims to offer an integrative and efficient way of performing RMSD calculations of large sets of structures. It is specially tuned to do fast collective RMSD calculations, as pairwise RMSD matrices.',
      author='Victor Alejandro Gil Sepulveda',
      author_email='victor.gil.sepulveda@gmail.com',
      url='https://github.com/victor-gil-sepulveda/pyRMSD.git',
      packages=['pyRMSD','pyRMSD.utils'],
      package_dir={'pyRMSD':'./pyRMSD'},
      py_modules=['pyRMSD.availableCalculators', 'pyRMSD.matrixHandler', 'pyRMSD.RMSDCalculator', 'pyRMSD.utils.proteinReading'],
      include_dirs = [numpy.get_include()],
      ext_modules=[
                   Extension('pyRMSD.pdbReader',[
                                          'src/pdbreaderlite/PDBReader.cpp',
                                          'src/python/readerLite.cpp'
                   ]),
                   Extension('pyRMSD.condensedMatrix', [
                                                 'src/matrix/Matrix.cpp',
                                                 'src/matrix/Statistics.cpp'
                   ]),
                   Extension(
                             'pyRMSD.calculators',
                             sources = [
                                        'src/python/pyRMSD.cpp',
                                        
                                        'src/calculators/RMSD.cpp',
                                        'src/calculators/RMSDTools.cpp',
                                        
                                        'src/calculators/QTRFIT/RMSDSerial.cpp',
                                        'src/calculators/QTRFIT/RMSDomp.cpp',

                                        'src/calculators/QCP/kernel_functions_omp.cpp',
                                        'src/calculators/QCP/kernel_functions_serial.cpp',
                                        'src/calculators/QCP/ThRMSDSerial.cpp',
                                        'src/calculators/QCP/ThRMSDSerialOmp.cpp'
                             ],
                             extra_compile_args=['-fopenmp'],
                             extra_link_args=['-fopenmp']
                   )
      ],
      
     )
