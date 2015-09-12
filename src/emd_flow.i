%module emd_flow
%{
#define SWIG_FILE_WITH_INIT
#include "python_helpers.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(const double* emd_costs, int num_emd_costs)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* data, int rows, int cols)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* support, int output_rows, int output_cols)}

%include "python_helpers.h"
