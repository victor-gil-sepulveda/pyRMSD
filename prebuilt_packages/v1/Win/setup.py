'''
Created on 25/02/2013

@author: victor
'''
from distutils.core import setup

setup(name='pyRMSD',
      version='1.0',
      description='pyRMSD is a small Python package that aims to offer an integrative and efficient way of performing RMSD calculations of large sets of structures. It is specially tuned to do fast collective RMSD calculations, as pairwise RMSD matrices.',
      author='Victor Alejandro Gil Sepulveda',
      author_email='victor.gil.sepulveda@gmail.com',
      url='https://github.com/victor-gil-sepulveda/pyRMSD.git',
      packages=['pyRMSD','pyRMSD.utils'],
      package_dir={'pyRMSD':'.'},
      package_data={'pyRMSD':['*.pyd']},
      py_modules=['availableCalculators', 'matrixHandler', 'RMSDCalculator', 'utils.proteinReading'],
      #data_files=[('.', ['calculators.pyd', 'condensedMatrix.pyd', 'pdbReader.pyd'])]
     )
