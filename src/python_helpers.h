#ifndef __EMDFLOW_PYTHON_HELPERS_H__
#define __EMDFLOW_PYTHON_HELPERS_H__

#include <vector>

#include "emd_flow_network_sap.h"

void emdflow_relaxed(const double* data, int rows, int cols,
                     const double* emd_costs, int num_emd_costs,
                     int sparsity,
                     double lambda,
                     int* support, int output_rows, int output_cols) {
  if (num_emd_costs <= 0) {
    throw std::invalid_argument("At least one EMD cost required.");
  }
  if (output_rows != rows || output_cols != cols) {
    throw std::invalid_argument("Output dimensions must match input "
        "dimensions.");
  }

  std::vector<std::vector<double> > input(rows);
  // numpy uses row major by default
  for (int ii = 0; ii < rows; ++ii) {
    input[ii].resize(cols);
    for (int jj = 0; jj < cols; ++jj) {
      input[ii][jj] = data[ii * cols + jj];
    }
  }
  std::vector<double> emd_costs2(num_emd_costs);
  for (int ii = 0; ii < num_emd_costs; ++ii) {
    emd_costs2[ii] = emd_costs[ii];
  }
  std::vector<std::vector<bool> > result;
  
  EMDFlowNetworkSAP algo(input, num_emd_costs - 1, emd_costs2);
  algo.set_sparsity(sparsity);
  algo.run_flow(lambda, 1.0);
  algo.get_support(&result);
  
  // numpy uses row major by default
  for (int ii = 0; ii < rows; ++ii) {
    for (int jj = 0; jj < cols; ++jj) {
      support[ii * cols + jj] = (result[ii][jj] ? 1 : 0);
    }
  }
}

#endif
