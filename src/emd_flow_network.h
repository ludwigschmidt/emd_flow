#ifndef __EMD_FLOW_NETWORK_H__
#define __EMD_FLOW_NETWORK_H__

#include <vector>
#include <string>

class EMDFlowNetwork {
 public:
  EMDFlowNetwork() { }
  virtual void set_sparsity(int k) = 0;
  virtual void run_flow(double lambda) = 0;
  virtual int get_EMD_used() = 0;
  virtual double get_supported_amplitude_sum() = 0;
  virtual void get_support(std::vector<std::vector<bool> >* support) = 0;
  virtual int get_num_nodes() = 0;
  virtual int get_num_edges() = 0;
  virtual int get_num_columns() = 0;
  virtual int get_num_rows() = 0;
  virtual void get_performance_diagnostics(std::string* s) { *s = "";}
  virtual ~EMDFlowNetwork() { }
};

#endif
