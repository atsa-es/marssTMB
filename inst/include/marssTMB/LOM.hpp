/// @file LOM.hpp
// Function for a list of matrices

#ifndef LOM_hpp
#define LOM_hpp 1

template <class Type>
struct LOM : vector<matrix<Type> > {
  LOM(SEXP x) {  // x = list passed from R
    (*this).resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

// List of integer vectors
template <class Type>
struct LOVi : vector<vector<int> > {
  LOVi(SEXP x) {  // x = list passed from R
    (*this).resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asVector<int>(sm);
    }
  }
};

// Convert the vectorized fixed, free, par in the non-vec matrix
template <class Type>
matrix<Type> parmat(matrix<Type> fixed, matrix<Type> free, 
            matrix<Type> par, vector<int> dims) {
  int nn = dims(0)*dims(1); /* num elements */
  int pp = par.col(0).size(); /* num parameters */
  matrix<Type> mat(dims(0), dims(1));
  matrix<Type> F(nn, pp);
  int i = 0;
  for (int r = 0; r < nn; r++)
    for (int c = 0; c < pp; c++)
      F(r,c) = free(i++,0);
  i = 0;
  vector<Type> fil = fixed + F*par; // vec of the matrix
  for (int r = 0; r < dims(1); r++)
    for (int c = 0; c < dims(0); c++)
      mat(r,c) = fil(i++);
  return mat;
};


#endif
