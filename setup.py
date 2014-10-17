"""
Created on 25/02/2013

@author: victor
"""
from distutils.core import setup, Extension
import numpy
import distutils.sysconfig
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
      name = 'pyRMSD',
      version = '4.1.0',
      description = 'pyRMSD is a small Python package that aims to offer an integrative and \
      efficient way of performing RMSD calculations of large sets of structures. It is specially \
      tuned to do fast collective RMSD calculations, as pairwise RMSD matrices.',
      author = 'Victor Alejandro Gil Sepulveda',
      author_email = 'victor.gil.sepulveda@gmail.com',
      url = 'https://github.com/victor-gil-sepulveda/pyRMSD.git',
      packages = [
                'pyRMSD',
                'pyRMSD.utils'
                ],
      package_dir = {'pyRMSD':'./pyRMSD'},
      py_modules = [
                  'pyRMSD.availableCalculators',
                  'pyRMSD.matrixHandler',
                  'pyRMSD.RMSDCalculator',
                  'pyRMSD.utils.proteinReading',
                  'pyRMSD.symmTools'
                  ],
      include_dirs = [
                      numpy.get_include(),
                      distutils.sysconfig.get_python_inc()
                      ],
      ext_modules = [
                   Extension('pyRMSD.pdbReader',[
                                          'src/pdbreaderlite/PDBReader.cpp',
                                          'src/pdbreaderlite/PDBReaderObject.cpp'
                   ]),
                   Extension('pyRMSD.condensedMatrix', [
                                                 'src/matrix/Matrix.cpp',
                                                 'src/matrix/Statistics.cpp'
                   ]),
                   Extension(
                             'pyRMSD.calculators',
                             sources = [
                                        'src/python/pyRMSD.cpp',

                                        'src/calculators/RMSDTools.cpp',
                                        'src/calculators/RMSDCalculator.cpp',
                                        'src/calculators/RMSDCalculationData.cpp',
                                        'src/calculators/KernelFunctions.cpp',

                                        'src/calculators/factory/RMSDCalculatorFactory.cpp',

                                        'src/calculators/KABSCH/KABSCHSerialKernel.cpp',
                                        'src/calculators/KABSCH/KABSCHOmpKernel.cpp',

                                        'src/calculators/QTRFIT/QTRFITSerialKernel.cpp',
                                        'src/calculators/QTRFIT/QTRFITOmpKernel.cpp',

                                        'src/calculators/QCP/QCPSerialKernel.cpp',
                                        'src/calculators/QCP/QCPSerialFloatKernel.cpp',
                                        'src/calculators/QCP/QCPOmpKernel.cpp',

                                        'src/calculators/NOSUP/NOSUPSerialKernel.cpp',
                                        'src/calculators/NOSUP/NOSUPOmpKernel.cpp',

                             ],
                             extra_compile_args = ['-fopenmp','-O3'],
                             extra_link_args = ['-fopenmp']
                   )
      ],
      license = 'LICENSE.txt',
      long_description = read('README.rst'),
      install_requires=[
        "numpy >= 1.6.1"
      ],
)
