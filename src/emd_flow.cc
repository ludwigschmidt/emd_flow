#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <memory>
#include <string>
#include <limits>

#include "emd_flow_network.h"
#include "emd_flow_network_factory.h"
#include "emd_flow.h"

using namespace std;

const int kOutputBufferSize = 10000;
char output_buffer[kOutputBufferSize];

// Make lambda larger until we find a solution that fits into the EMD budget.
// Returns true if we find a solution in [emd_bound_low, emd_bound_high].
bool increase_lambda(const emd_flow_args& args, emd_flow_result* result,
    EMDFlowNetwork* network, double* lambda_high);

// Make lambda smaller until we find a solution that does not fit into the
// EMD budget. Return true if we find a solution in
// [emd_bound_low, emd_bound_high].
bool decrease_lambda(const emd_flow_args& args, double current_lambda_high,
    emd_flow_result* result, EMDFlowNetwork* network, double* lambda_low);

// Binary search over lambda
void binary_search_lambda(const emd_flow_args& args, double lambda_low,
    double lambda_high, emd_flow_result* result, EMDFlowNetwork* network);

// Set the result struct to values indicating an error.
void clear_result(emd_flow_result* result);


// Main routine for computing the EMD flow.
void emd_flow(const emd_flow_args& args, emd_flow_result* result) {
  clock_t total_time_begin = clock();

  int r = args.x.size();
  int c = args.x[0].size();

  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "r = %d,  c = %d,  s = %d,  "
        "emd_bound_low = %d, emd_bound_high = %d\n", r, c, args.s,
        args.emd_bound_low, args.emd_bound_high);
    args.output_function(output_buffer);
    snprintf(output_buffer, kOutputBufferSize, "lambda_low = %e, "
        "lambda_high = %e, num_search_iterations = %d\n", args.lambda_low,
        args.lambda_high, args.num_search_iterations);
    args.output_function(output_buffer);
  }

  /*for (size_t ii = 0; ii < r; ++ii) {
    for (size_t jj = 0; jj < c; ++jj) {
      printf("%lf ", args.x[ii][jj]);
    }
    printf("\n");
  }
  printf("\n");*/

  // build graph
  clock_t graph_construction_time_begin = clock();

  int outdegree_vertical_distance = args.outdegree_vertical_distance;
  if (outdegree_vertical_distance == -1) {
    outdegree_vertical_distance = r - 1;
  } else if (outdegree_vertical_distance < -1) {
    snprintf(output_buffer, kOutputBufferSize, "Error: "
        "outdegree_vertical_distance cannot be less than -1, given value is %d"
        ".\n", outdegree_vertical_distance);
    args.output_function(output_buffer);
    clear_result(result);
    return;
  }

  vector<double> emd_costs = args.emd_costs;
  if (emd_costs.size() == 0) {
    for (int ii = 0; ii <= outdegree_vertical_distance; ++ii) {
      emd_costs.push_back(ii);
    }
  } else if (static_cast<int>(emd_costs.size())
      != outdegree_vertical_distance + 1) {
    snprintf(output_buffer, kOutputBufferSize, "Error: "
        "the emd_costs vector has an incorrect number of entries: %lu entries, "
        " should be %d\n", emd_costs.size(), outdegree_vertical_distance + 1);
    args.output_function(output_buffer);
    clear_result(result);
    return;
  }

  auto_ptr<EMDFlowNetwork> network =
      EMDFlowNetworkFactory::create_EMD_flow_network(args.x,
      outdegree_vertical_distance, emd_costs, args.alg_type);
  network->set_sparsity(args.s);

  clock_t graph_construction_time = clock() - graph_construction_time_begin;

  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "The graph has %d nodes and %d "
        "edges.\n", network->get_num_nodes(), network->get_num_edges());
    args.output_function(output_buffer);
    snprintf(output_buffer, kOutputBufferSize, "Total construction time: %f "
        "s\n ", static_cast<double>(graph_construction_time) / CLOCKS_PER_SEC);
    args.output_function(output_buffer);
  }

  double lambda_high = args.lambda_high;
  double lambda_low = args.lambda_low;

  if (!increase_lambda(args, result, network.get(), &lambda_high)) {
    if (!decrease_lambda(args, lambda_high, result, network.get(),
        &lambda_low)) {
      binary_search_lambda(args, lambda_low, lambda_high, result,
          network.get());
    }
  }

  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Final l: %e, amp sum: %e, "
        "EMD cost: %d\n", lambda_high, result->amp_sum, result->emd_cost);
    args.output_function(output_buffer);
  }

  clock_t total_time = clock() - total_time_begin;
  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Total time %f s\n",
        static_cast<double>(total_time) / CLOCKS_PER_SEC);
    args.output_function(output_buffer);

    string performance_diagnostics;
    network->get_performance_diagnostics(&performance_diagnostics);
    snprintf(output_buffer, kOutputBufferSize, "Performance diagnostics:\n%s\n",
        performance_diagnostics.c_str());
    args.output_function(output_buffer);
  }

  return;
}

void binary_search_lambda(const emd_flow_args& args, double lambda_low,
    double lambda_high, emd_flow_result* result, EMDFlowNetwork* network) {

  // binary search on lambda
  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Binary search on lambda ...\n");
    args.output_function(output_buffer);
  }

  int cur_emd_cost = 0;
  double cur_amp_sum = 0;
  int current_iteration = 1;

  while (current_iteration <= args.num_search_iterations
      && (cur_emd_cost < args.emd_bound_low
      || cur_emd_cost > args.emd_bound_high)) {
    ++current_iteration;
    double cur_lambda = (lambda_high + lambda_low) / 2;
    network->run_flow(cur_lambda, 1.0);
    cur_emd_cost = network->get_EMD_used();
    cur_amp_sum = network->get_supported_amplitude_sum();

    if (args.verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l_cur: %e  (l_low: %e, "
          "l_high: %e)  EMD: %d  amp sum: %e\n", cur_lambda, lambda_low,
          lambda_high, cur_emd_cost, cur_amp_sum);
      args.output_function(output_buffer);
    }

    if (cur_emd_cost <= args.emd_bound_high) {
      lambda_high = cur_lambda;
    } else {
      lambda_low = cur_lambda;
    }
  }

  // TODO: don't rerun flow here?
  // run with final lambda
  network->run_flow(lambda_high, 1.0);
  result->final_lambda_low = lambda_low;
  result->final_lambda_high = lambda_high;
  result->emd_cost = network->get_EMD_used();
  result->amp_sum = network->get_supported_amplitude_sum();
  network->get_support(result->support);
}

// Increase lambda until we find a lambda such that the EMD cost is smaller
// than the lower EMD bound. If we find a feasible solution during the process,
// we return true. Otherwise we return false.
bool increase_lambda(const emd_flow_args& args, emd_flow_result* result,
    EMDFlowNetwork* network, double* lambda_high) {
  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize,
        "Finding large enough value of lambda ...\n");
    args.output_function(output_buffer);
  }

  // Check what the cheapest flow is (ignoring node costs). If even the
  // cheapest flow costs more than the upper EMD bound, we cannot satisfy it.
  network->run_flow(1.0, 0.0);
  int cur_emd_cost = network->get_EMD_used();
  double cur_amp_sum = network->get_supported_amplitude_sum();

  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "l_EMD: 1.0  l_signal: 0.0  "
        "EMD: %d  amp sum: %e\n", cur_emd_cost, cur_amp_sum);
    args.output_function(output_buffer);
  }

  if (cur_emd_cost > args.emd_bound_high) {
    snprintf(output_buffer, kOutputBufferSize, "Cannot satisfy the upper EMD "
        "bound regardless of the signal approximation: the smallest feasible "
        "EMD cost is %d while the upper EMD bound is %d. Consider changing the "
        "EMD bounds or the edge EMD costs.", cur_emd_cost, args.emd_bound_high);
    args.output_function(output_buffer);
    clear_result(result);
    return true;
  }

  while (true) {
    network->run_flow(*lambda_high, 1.0);
    cur_emd_cost = network->get_EMD_used();
    cur_amp_sum = network->get_supported_amplitude_sum();

    if (args.verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l: %e  EMD: %d  amp sum: %e"
          "\n", *lambda_high, cur_emd_cost, cur_amp_sum);
      args.output_function(output_buffer);
    }

    if (cur_emd_cost <= args.emd_bound_high) {
      if (cur_emd_cost >= args.emd_bound_low) {
        result->final_lambda_low = args.lambda_low;
        result->final_lambda_high = *lambda_high;
        result->emd_cost = cur_emd_cost;
        result->amp_sum = cur_amp_sum;
        network->get_support(result->support);
        return true;
      } else {
        return false;
      }
    } else {
      *lambda_high = *lambda_high * 2;
    }
  }
}

// Decrease lambda until we find a lambda such that the EMD cost is larger
// than the upper EMD bound. If we find a feasible solution during the process,
// we return true. Otherwise we return false.
bool decrease_lambda(const emd_flow_args& args, double current_lambda_high,
    emd_flow_result* result, EMDFlowNetwork* network, double* lambda_low) {
  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize,
        "Finding small enough value of lambda ...\n");
    args.output_function(output_buffer);
  }

  // Calculate the best approximation that satisfies only s-sparsity in
  // each column. This allows us to end early in case the EMD bounds are
  // larger than what we need for a perfect approximation.
  *lambda_low = 0.0;
  network->run_flow(*lambda_low, 1.0);
  int cur_emd_cost = network->get_EMD_used();
  double cur_amp_sum = network->get_supported_amplitude_sum();
  if (args.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "l: %e  EMD: %d  amp sum: %e"
        "\n", *lambda_low, cur_emd_cost, cur_amp_sum);
    args.output_function(output_buffer);
  }

  // In this case, the final EMD cost might be less than emd_bound_low.
  // But we are running with lambda=0, so we cannot use more EMD.
  if (cur_emd_cost < args.emd_bound_high) {
    if (args.emd_bound_low < args.emd_bound_high
        && cur_emd_cost < args.emd_bound_low) {
      snprintf(output_buffer, kOutputBufferSize, "Found a solution with lambda "
          "= 0, so the solution does not satisfy the lower EMD bound.");
      args.output_function(output_buffer);
    }

    result->final_lambda_low = *lambda_low;
    result->final_lambda_high = current_lambda_high;
    result->emd_cost = cur_emd_cost;
    result->amp_sum = cur_amp_sum;
    network->get_support(result->support);
    return true;
  }

  *lambda_low = args.lambda_low;
  while (true) {
    network->run_flow(*lambda_low, 1.0);
    cur_emd_cost = network->get_EMD_used();
    cur_amp_sum = network->get_supported_amplitude_sum();

    if (args.verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l: %e  EMD: %d  amp sum: %e"
          "\n", *lambda_low, cur_emd_cost, cur_amp_sum);
      args.output_function(output_buffer);
    }

    if (cur_emd_cost > args.emd_bound_high) {
      break;
    } else {
      if (cur_emd_cost >= args.emd_bound_low) {
        result->final_lambda_low = *lambda_low;
        result->final_lambda_high = current_lambda_high;
        result->emd_cost = cur_emd_cost;
        result->amp_sum = cur_amp_sum;
        network->get_support(result->support);
        return true;
      }
      *lambda_low = *lambda_low / 2;
    }
  }
  return false;
}

void clear_result(emd_flow_result* result) {
  result->support->clear();
  result->emd_cost = 0;
  result->amp_sum = 0;
  result->final_lambda_low = 0;
  result->final_lambda_high = 0;
}
