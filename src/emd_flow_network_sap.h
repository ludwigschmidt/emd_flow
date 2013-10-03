#ifndef __EMD_FLOW_NETWORK_SAP_H__
#define __EMD_FLOW_NETWORK_SAP_H__

#include "emd_flow_network.h"

#include <algorithm>
#include <vector>
#include <cstddef>

class EMDFlowNetworkSAP : public EMDFlowNetwork {
 public:
  EMDFlowNetworkSAP(
      const std::vector<std::vector<double> >& amplitudes,
      int outdegree_vertical_distance,
      const std::vector<double>& emd_costs);
  void set_sparsity(int s);
  void run_flow(double EMD_lambda, double signal_lambda);
  int get_EMD_used();
  double get_supported_amplitude_sum();
  void get_support(std::vector<std::vector<bool> >* support);
  int get_num_nodes();
  int get_num_edges();
  int get_num_columns();
  int get_num_rows();
  void get_performance_diagnostics(std::string* s);
  ~EMDFlowNetworkSAP() { }

 private:
  // node indices:
  // source: 0
  // sink: 1
  // other innode: 2 + 2 * (c_ * num_rows + r_)
  // other outnode: 3 + 2 * (c_ * num_rows + r_)
  typedef size_t NodeIndex;
  typedef size_t EdgeIndex;

  struct Edge {
    NodeIndex to;
    int capacity;
    double cost;
    EdgeIndex opposite;

    Edge(NodeIndex _to, int _capacity, double _cost, EdgeIndex _opposite)
        : to(_to), capacity(_capacity), cost(_cost), opposite(_opposite) { }
  };

  // amplitudes
  std::vector<std::vector<double> > a_;
  // sparsity per column
  int sparsity_;
  // number of rows
  int r_;
  // number of columns
  int c_;
  // maximum vertical distance covered by an edge between neighboring columns
  int outdegree_vertical_distance_;
  // emd costs for an edge between columns with vertical distance i.
  std::vector<double> emd_costs_;

  // source, sink
  NodeIndex s_, t_;

  // edges representing a node cost
  std::vector<std::vector<EdgeIndex> > node_edges_;
  // edges representing an EMD step
  std::vector<std::vector<std::vector<EdgeIndex> > > emd_edges_;
  // edges leaving a node
  std::vector<std::vector<EdgeIndex> > outgoing_edges_;
  // set of all edges
  std::vector<Edge> e_;

  // node potentials
  std::vector<double> potential_;

  long long total_inner_iterations;
  long long checking_inner_iterations;
  long long updating_inner_iterations;

  NodeIndex innode_index(int r, int c) {
    return 2 + 2 * (c * r_ + r);
  }

  NodeIndex outnode_index(int r, int c) {
    return innode_index(r, c) + 1;
  }

  size_t num_destinations(int r) {
    return 1 + std::min(outdegree_vertical_distance_, r)
             + std::min(outdegree_vertical_distance_, r_ - r - 1);
  }

  int first_destination(int r) {
    return std::max(0, r - outdegree_vertical_distance_);
  }

  void apply_EMD_lambda(double lambda);
  void apply_signal_lambda(double lambda);
  void reset_flow();
  void compute_initial_potential();
  void print_full_graph();
};


#endif
