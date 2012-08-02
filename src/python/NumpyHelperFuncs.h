/////////////////////////////////////////////////////////////////////////////
/// \file NumpyHelperFuncs.h
/// \brief bla, bla
///
/// \author vgil
/// \date 22/02/2012
/////////////////////////////////////////////////////////////////////////////

#ifndef NUMPYHELPERFUNCS_H_
#define NUMPYHELPERFUNCS_H_

class PyArrayObject;

int      not_doublevector(PyArrayObject *vec);
double** ptrvector(long n);
double*  pyvector_to_Carrayptrs(PyArrayObject *arrayin);

#endif /* NUMPYHELPERFUNCS_H_ */
