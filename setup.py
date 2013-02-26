'''
Created on 25/02/2013

@author: victor
'''
from distutils.core import setup, Extension

setup(name='pyRMSD',
      version='1.0',
      description='pyRMSD is a small Python package that aims to offer an integrative and efficient way of performing RMSD calculations of large sets of structures. It is specially tuned to do fast collective RMSD calculations, as pairwise RMSD matrices.',
      author='Victor Alejandro Gil Sepulveda',
      author_email='victor.gil.sepulveda@gmail.com',
      url='https://github.com/victor-gil-sepulveda/pyRMSD.git',
      packages=['pyRMSD','pyRMSD.utils'],
      package_dir={'pyRMSD':'../pyRMSD'},
      py_modules=['availableCalculators', 'matrixHandler', 'RMSDCalculator', 'utils.proteinReading'],
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
                                        
                                        'src/serial/RMSD.cpp',
                                        'src/serial/RMSDomp.cpp',
                                        'src/serial/RMSDSerial.cpp',
                                        'src/serial/RMSDTools.cpp',
                                        
                                        'src/theobald/kernel_functions_omp.cpp',
                                        'src/theobald/kernel_functions_serial.cpp',
                                        'src/theobald/ThRMSDSerial.cpp',
                                        'src/theobald/ThRMSDSerialOmp.cpp'
                             ],
                             extra_compile_args=['-fopenmp'],
                             extra_link_args=['-fopenmp']
                   )
      ],
      
     )
