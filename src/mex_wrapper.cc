#include <vector>
#include <string>
#include <set>

#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "mex_helper.h"
#include "emd_flow.h"
#include "emd_flow_network_factory.h"

using namespace std;

void output_function(const char* s) {
  mexPrintf(s);
  mexEvalString("drawnow;");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 3) {
    mexErrMsgTxt("At least three input argument required (amplitudes, sparsity,"
        " EMD budget.");
  }
  if (nrhs > 4) {
    mexErrMsgTxt("Too many input arguments, at most four: amplitudes, sparsity,"
        " EMD budget and the options struct.");
  }
  if (nlhs > 5) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  vector<vector<double> > a;
  if (!get_double_matrix(prhs[0], &a)) {
    mexErrMsgTxt("Amplitudes need to be a two-dimensional double array.");
  }

  int s = 0;
  if (!get_double_as_int(prhs[1], &s)) {
    mexErrMsgTxt("Sparsity has to be a double scalar.");
  }

  int emd_bound_low = 0;
  int emd_bound_high = 0;
  if (!get_double_as_int(prhs[2], &emd_bound_low)) {
    if (!get_double_interval_as_ints(prhs[2], &emd_bound_low,
        &emd_bound_high)) {
      mexErrMsgTxt("EMD budget has to be a double scalar or a double "
          "interval.");
    }
  } else {
    emd_bound_high = emd_bound_low;
  }

  // optional parameters
  bool verbose = false;
  double lambda_low = 0.5;
  double lambda_high = 1.0;
  int num_iter = 10;
  vector<double> emd_costs;
  if (nrhs == 4) {
    set<string> known_options;
    known_options.insert("verbose");
    known_options.insert("lambda_low");
    known_options.insert("lambda_high");
    known_options.insert("num_iterations");
    vector<string> options;
    if (!get_fields(prhs[3], &options)) {
      mexErrMsgTxt("Cannot get fields from options argument.");
    }
    for (size_t ii = 0; ii < options.size(); ++ii) {
      if (known_options.find(options[ii]) == known_options.end()) {
        const size_t tmp_size = 1000;
        char tmp[tmp_size];
        snprintf(tmp, tmp_size, "Unknown option \"%s\"\n", options[ii].c_str());
        mexErrMsgTxt(tmp);
      }
    }

    if (has_field(prhs[3], "verbose")
        && !get_bool_field(prhs[3], "verbose", &verbose)) {
      mexErrMsgTxt("verbose flag has to be a boolean scalar.");
    }

    if (has_field(prhs[3], "lambda_low")
        && !get_double_field(prhs[3], "lambda_low", &lambda_low)) {
      mexErrMsgTxt("lambda_low field has to be a double scalar.");
    }

    if (has_field(prhs[3], "lambda_high")
        && !get_double_field(prhs[3], "lambda_high", &lambda_high)) {
      mexErrMsgTxt("lambda_high field has to be a double scalar.");
    }

    if (has_field(prhs[3], "num_iterations")
        && !get_double_field_as_int(prhs[3], "num_iterations", &num_iter)) {
      mexErrMsgTxt("num_iterations field has to be a double scalar.");
    }
  }

  emd_flow_args args(a);
  args.s = s;
  args.emd_bound_low = emd_bound_low;
  args.emd_bound_high = emd_bound_high;
  args.lambda_low = lambda_low;
  args.lambda_high = lambda_high;
  args.num_search_iterations = num_iter;
  args.outdegree_vertical_distance = -1;
  args.emd_costs = emd_costs;
  args.alg_type = EMDFlowNetworkFactory::kShortestAugmentingPath;
  args.output_function = output_function;
  args.verbose = verbose;

  
  emd_flow_result result;

  emd_flow(args, &result);

  if (nlhs >= 1) {
    set_double_matrix(&(plhs[0]), *(result.support));
  }

  if (nlhs >= 2) {
    set_double(&(plhs[1]), result.emd_cost);
  }

  if (nlhs >= 3) {
    set_double(&(plhs[2]), result.amp_sum);
  }

  if (nlhs >= 4) {
    set_double(&(plhs[3]), result.final_lambda_low);
  }

  if (nlhs >= 5) {
    set_double(&(plhs[4]), result.final_lambda_high);
  }

  /*
  mexPrintf("r = %d, c = %d, k = %d, EMD budget = %d\n", r, c, k, emd_budget);
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < c; ++jj) {
      mexPrintf("%lf ", a[ii][jj]);
    }
    mexPrintf("\n");
  }
  */
}
