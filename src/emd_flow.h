#ifndef __EMD_FLOW_H__
#define __EMD_FLOW_H__

#include <vector>

#include "emd_flow_network_factory.h"

struct emd_flow_args {
  // input amplitudes (will not be squared)
  const std::vector<std::vector<double> >& x;
  // sparsity per column
  int s;
  // Bounds on the EMD budget. If we find a solution with EMD cost
  // between the lower bound and the upper bound we return it.
  int emd_bound_low;
  int emd_bound_high;
  // Initial guesses for bounds on lambda
  // The guesses are corrected if the real lambda is not between them.
  double lambda_low;
  double lambda_high;
  // Maximum number of search iterations
  int num_search_iterations;
  // Maximum vertical distance covered by edges from one columns to the next.
  // -1 indicates that we should build a full graph.
  int outdegree_vertical_distance;
  // Cost of the edges from one column to the next.
  // The vector should have outdegree_vertical_distance + 1 (or num_rows, if
  // outdegree_vertical_distance is -1) entries. Entry i (starting at 0) should
  // contain the cost of an edge with vertical distance i from one column to
  // the next. Setting the vector to [0, 1, 2, 3, ..., ] gives the standard EMD.
  // Passing an empty vector defaults to the EMD.
  std::vector<double> emd_costs;
  // The internal flow algorithm to use
  EMDFlowNetworkFactory::EMDFlowNetworkType alg_type;
  // The output function
  void (*output_function)(const char*);
  // Verbose output?
  bool verbose; 

  emd_flow_args(const std::vector<std::vector<double> >& x_): x(x_) { }
};

struct emd_flow_result {
  // Bool matrix indicating the support (same dimensions as input)
  std::vector<std::vector<bool> >* support;
  // EMD cost of the solution
  int emd_cost;
  // Absolute amplitude sum of the solution
  double amp_sum;
  // Final values of bounds on lambda
  double final_lambda_low;
  double final_lambda_high;
};

void emd_flow(
    const emd_flow_args& args,
    emd_flow_result* result);

#endif
