#include "emd_flow.h"
#include "emd_flow_network.h"

#include <cstdio>                                                               
#include <vector>

#include "boost/assign/list_of.hpp"
#include "gtest/gtest.h"

using namespace boost::assign;
using namespace std;

void WriteToStderr(const char* s) {
  fprintf(stderr, s);
  fflush(stderr);
}

void FillArgs(int s, int B, emd_flow_args* args) {
  args->s = s;
  args->emd_bound_low = B;
  args->emd_bound_high = B;
  args->lambda_low = 0.5;
  args->lambda_high = 1;
  args->num_search_iterations = 10;
  args->outdegree_vertical_distance = -1;
  args->alg_type = EMDFlowNetworkFactory::kShortestAugmentingPath;
  args->output_function = WriteToStderr;
  args->verbose = true;
}

void CheckResultConsistency(const emd_flow_result&) {
  // TODO:implement
}

void CheckResult(const emd_flow_result& result,
                 int expected_EMD,
                 double expected_amplitude_sum) {
  EXPECT_EQ(result.emd_cost, expected_EMD);
  EXPECT_DOUBLE_EQ(result.amp_sum, expected_amplitude_sum);
  EXPECT_LT(result.final_lambda_low, result.final_lambda_high);
}

void CheckResult(const emd_flow_result& result,
                 const vector<vector<bool > > expected_support,
                 int expected_EMD,
                 double expected_amplitude_sum) {
  ASSERT_EQ(result.support->size(), expected_support.size());
  for (size_t ii = 0; ii < expected_support.size(); ++ii) {
    ASSERT_EQ((*result.support)[ii].size(), expected_support[ii].size());
    for (size_t jj = 0; jj < expected_support[ii].size(); ++jj) {
      EXPECT_EQ((*result.support)[ii][jj], expected_support[ii][jj])
          << "support mismatch at (" << ii << ", " << jj << ")";
    }
  }
  CheckResult(result, expected_EMD, expected_amplitude_sum);
}

TEST(EMDFlowTest, SimpleZeroEMDTwoSparsity) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(1.0));
  x.push_back(list_of(1.0)(0.0));
  x.push_back(list_of(0.0)(0.0));
  const int s = 2;
  const int B = 0;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(1)(1));
  expected_support.push_back(list_of(1)(1));
  expected_support.push_back(list_of(0)(0));
  CheckResult(result, expected_support, 0, 2.0);
}

TEST(EMDFlowTest, SimpleZeroEMDOneSparsity) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(101.0));
  x.push_back(list_of(100.0)(0.0));
  x.push_back(list_of(0.0)(0.0));
  const int s = 1;
  const int B = 0;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(1)(1));
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(0)(0));
  CheckResult(result, expected_support, 0, 101.0);
}

TEST(EMDFlowTest, SimpleTwoEMDOneSparsity) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(100.0));
  x.push_back(list_of(0.0)(0.0));
  x.push_back(list_of(100.0)(0.0));
  const int s = 1;
  const int B = 2;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(0)(1));
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(1)(0));
  CheckResult(result, expected_support, 2, 200.0);
}

TEST(EMDFlowTest, SimpleOneEMDOneSparsity) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(100.0));
  x.push_back(list_of(0.0)(0.0));
  x.push_back(list_of(101.0)(0.0));
  const int s = 1;
  const int B = 1;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(1)(1));
  CheckResult(result, expected_support, 0, 101.0);
}

TEST(EMDFlowTest, SimpleOneEMDOneSparsity2) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(1.1));
  x.push_back(list_of(0.0)(0.0));
  x.push_back(list_of(1.0)(0.0));
  const int s = 1;
  const int B = 1;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(1)(1));
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(0)(0));
  CheckResult(result, expected_support, 0, 1.1);
}

TEST(EMDFlowTest, SimpleOneEMDOneSparsity3) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(1.1));
  x.push_back(list_of(0.0)(1.0));
  x.push_back(list_of(1.0)(0.0));
  const int s = 1;
  const int B = 1;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(0)(1));
  expected_support.push_back(list_of(1)(0));
  CheckResult(result, expected_support, 1, 2.0);
}

TEST(EMDFlowTest, SimpleOneEMDOneSparsity4) {
  vector<vector<double> > x;
  x.push_back(list_of(0.0)(1.0));
  x.push_back(list_of(0.0)(0.25));
  x.push_back(list_of(10.0)(0.0));
  const int s = 1;
  const int B = 1;
  emd_flow_args args(x);
  FillArgs(s, B, &args);

  vector<vector<bool> > support;
  emd_flow_result result;
  result.support = &support;
  emd_flow(args, &result); 

  vector<vector<bool> > expected_support;
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(0)(0));
  expected_support.push_back(list_of(1)(1));
  CheckResult(result, expected_support, 0, 10.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);                                       
  return RUN_ALL_TESTS();
}
