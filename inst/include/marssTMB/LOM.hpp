/// @file LOM.hpp
// Function for a list of matrices

#ifndef LOM_hpp
#define LOM_hpp 1

// List of Matrices
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

// Convert the vectorized fixed, free, par into a 2D matrix
// fixed, free, and par are all column vectors, dims is a vector of the
// matrix dimensions
template <class Type>
matrix<Type> parmat2(matrix<Type> fixed, matrix<Type> free, 
                    matrix<Type> par, vector<int> dims) {
  int vecdim = dims(0)*dims(1);
  matrix<Type> F(vecdim, par.rows());
  free.resize(vecdim, par.rows());
  matrix<Type> mat(vecdim, 1);
  mat = fixed + free*par;
  mat.resize(dims(0), dims(1));
  return mat;
};

// Convert the vectorized fixed, free, par into matrix with each column
// a vec of free * par and time along the columns
template <class Type>
matrix<Type> parvec(matrix<Type> fixed, matrix<Type> free, 
                    matrix<Type> par, vector<int> dims, 
                    int tfree, int tfixed, int timeSteps) {
  matrix<Type> mat(dims(0)*dims(1), timeSteps);
  int rowD = dims(0)*dims(1);
  matrix<Type> rowOne(1, timeSteps); 
  rowOne.setOnes();
  if(tfixed == 1) fixed = fixed * rowOne;
  if(tfree == 1) free = free * rowOne;
  matrix<Type> I(rowD,rowD);
  I.setIdentity();
  matrix<Type> tpar(1,par.rows());
  tpar = par.transpose();
  mat = fixed + tmbutils::kronecker(tpar, I) * free;
  return mat;
};

// Gets the parameters associated with elem and returns a col vector (n x 1)
template <class Type>
matrix<Type> par(vector<Type> pars, vector<int> numpar, int elem) {
  if(numpar(elem) == 0){
    matrix<Type> mat(1, 1);
    mat.setZero();
    return mat;
  }else{
    matrix<Type> mat(numpar(elem), 1);
    mat.setZero();
    int i = 0;
    if(elem > 0) for (int r = 0; r < elem; r++) i = i + numpar(r);
    for (int z = i; z < (i + numpar(elem)); z++) mat(z-i,0) = pars(z);
    return mat;
  }
};


#endif
