/// @file isNA.hpp
// Function for detecting NAs

#ifndef isNA_hpp
#define isNA_hpp 1

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

#endif
