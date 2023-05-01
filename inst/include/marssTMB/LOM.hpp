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

#endif
