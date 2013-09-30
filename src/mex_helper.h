#ifndef __MEX_HELPER_H__
#define __MEX_HELPER_H__

#include <mex.h>
#include <matrix.h>
#include <cmath>
#include <vector>
#include <string>

bool get_double(const mxArray* raw_data, double* data) {
  int numdims = mxGetNumberOfDimensions(raw_data);
  const mwSize* dims = mxGetDimensions(raw_data);
  if (numdims != 2 || dims[0] != 1 || dims[1] != 1) {
    return false;
  }
  if (!mxIsClass(raw_data, "double")) {
    return false;
  }
  *data = static_cast<double*>(mxGetData(raw_data))[0];
  return true;
}

bool get_double_as_int(const mxArray* raw_data, int* data) {
  double tmp;
  if (get_double(raw_data, &tmp)) {
    *data = static_cast<double>(round(tmp) + 0.1);
    return true;
  } else {
    return false;
  }
}

bool get_bool(const mxArray* raw_data, bool* data) {
  int numdims = mxGetNumberOfDimensions(raw_data);
  const mwSize* dims = mxGetDimensions(raw_data);
  if (numdims != 2 || dims[0] != 1 || dims[1] != 1) {
    return false;
  }
  if (!mxIsClass(raw_data, "logical")) {
    return false;
  }
  *data = static_cast<bool*>(mxGetData(raw_data))[0];
  return true;
}

bool get_double_row_vector(const mxArray* raw_data,
    std::vector<double>* data) {
  int numdims = mxGetNumberOfDimensions(raw_data);
  const mwSize* dims = mxGetDimensions(raw_data);
  if (numdims != 2) {
    return false;
  }
  if (!mxIsClass(raw_data, "double")) {
    return false;
  }
  size_t r = dims[0];
  size_t c = dims[1];
  if (r != 1) {
    return false;
  }
  double* data_linear = static_cast<double*>(mxGetData(raw_data));
  if (data->size() != c) {
    data->resize(c);
  }
  for (size_t ic = 0; ic < c; ++ic) {
    (*data)[ic] = data_linear[ic];
  }
  return true;
}

bool get_double_interval_as_ints(const mxArray* raw_data, int* left,
    int* right) {
  std::vector<double> tmp;
  get_double_row_vector(raw_data, &tmp);
  if (tmp.size() != 2) {
    return false;
  }
  *left = static_cast<double>(round(tmp[0]) + 0.1);
  *right= static_cast<double>(round(tmp[1]) + 0.1);
  return true;
}

bool get_double_matrix(const mxArray* raw_data,
    std::vector<std::vector<double> >* data) {
  int numdims = mxGetNumberOfDimensions(raw_data);
  const mwSize* dims = mxGetDimensions(raw_data);
  if (numdims != 2) {
    return false;
  }
  if (!mxIsClass(raw_data, "double")) {
    return false;
  }
  size_t r = dims[0];
  size_t c = dims[1];
  double* data_linear = static_cast<double*>(mxGetData(raw_data));
  if (data->size() != r) {
    data->resize(r);
  }
  for (size_t ir = 0; ir < r; ++ir) {
    if ((*data)[ir].size() != c) {
      (*data)[ir].resize(c);
    }
    for (size_t ic = 0; ic < c; ++ic) {
      (*data)[ir][ic] = data_linear[ir + ic * r];
    }
  }
  return true;
}

bool get_fields(const mxArray* struc, std::vector<std::string>* fields) {
  if (!mxIsStruct(struc)) {
    return false;
  }
  size_t num_elems = mxGetNumberOfFields(struc);
  if (fields->size() != num_elems) {
    fields->resize(num_elems);
  }
  for (size_t ii = 0; ii < num_elems; ++ii) {
    (*fields)[ii] = mxGetFieldNameByNumber(struc, ii);
  }
  return true;
}

bool has_field(const mxArray* struc, const char* name) {
  if (!mxIsStruct(struc)) {
    return false;
  }
  mxArray* raw_data = mxGetField(struc, 0, name);
  return (raw_data != NULL);
}

bool get_double_field(const mxArray* struc, const char* name, double* data) {
  if (!mxIsStruct(struc)) {
    return false;
  }
  mxArray* raw_data = mxGetField(struc, 0, name);
  if (raw_data == NULL) {
    return false;
  }
  return get_double(raw_data, data);
}

bool get_double_field_as_int(const mxArray* struc, const char* name,
    int* data) {
  if (!mxIsStruct(struc)) {
    return false;
  }
  mxArray* raw_data = mxGetField(struc, 0, name);
  if (raw_data == NULL) {
    return false;
  }
  return get_double_as_int(raw_data, data);
}

bool get_bool_field(const mxArray* struc, const char* name, bool* data) {
  if (!mxIsStruct(struc)) {
    return false;
  }
  mxArray* raw_data = mxGetField(struc, 0, name);
  if (raw_data == NULL) {
    return false;
  }
  return get_bool(raw_data, data);
}

void set_double(mxArray** raw_data, double data) {
  *raw_data = mxCreateDoubleMatrix(1, 1, mxREAL);
  *(static_cast<double*>(mxGetData(*raw_data))) = data;
}

void set_double_matrix(mxArray** raw_data,
    const std::vector<std::vector<double> >& data) {
  int numdims = 2;
  mwSize dims[2];
  size_t r = data.size();
  size_t c = data[0].size();
  dims[0] = r;
  dims[1] = c;
  *raw_data = mxCreateNumericArray(numdims, dims, mxDOUBLE_CLASS, mxREAL);
  double* result_linear = static_cast<double*>(mxGetData(*raw_data));

  for (size_t ir = 0; ir < r; ++ir) {
    for (size_t ic = 0; ic < c; ++ic) {
      result_linear[ir + ic * r] = data[ir][ic];
    }
  }
}
void set_double_matrix(mxArray** raw_data,
    const std::vector<std::vector<bool> >& data) {
  std::vector<std::vector<double> > tmp_data;
  size_t r = data.size();
  size_t c = data[0].size();
  tmp_data.resize(r);
  for (size_t ir = 0; ir < r; ++ir) {
    tmp_data[ir].resize(c);
    for (size_t ic = 0; ic < c; ++ic) {
      tmp_data[ir][ic] = data[ir][ic];
    }
  }
  set_double_matrix(raw_data, tmp_data);
}

#endif
