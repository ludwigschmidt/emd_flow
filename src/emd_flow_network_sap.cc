#include "emd_flow_network_sap.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <queue>
#include <map>

using namespace std;

EMDFlowNetworkSAP::EMDFlowNetworkSAP(
    const std::vector<std::vector<double> >& amplitudes,
    int outdegree_vertical_distance,
    const vector<double>& emd_costs)
      : a_(amplitudes),
        outdegree_vertical_distance_(outdegree_vertical_distance),
        emd_costs_(emd_costs),
        total_inner_iterations(0),
        checking_inner_iterations(0),
        updating_inner_iterations(0) {
  r_ = amplitudes.size();
  c_ = amplitudes[0].size();

  a_.resize(r_);
  for (int row = 0; row < r_; ++row) {
    a_[row].resize(c_);
    for (int col = 0; col < c_; ++col) {
      a_[row][col] = amplitudes[row][col];
    }
  }

  // source and sink
  s_ = 0;
  t_ = 1;

  // potentials
  int num_nodes = 2 + 2 * r_ * c_;
  outgoing_edges_.resize(num_nodes);
  potential_.resize(outgoing_edges_.size());

  // add arcs from source to column 1
  for (int ii = 0; ii < r_; ++ii) {
    EdgeIndex next_edge_index = e_.size();
    Edge forward(innode_index(ii, 0), 1, 0.0, next_edge_index + 1);
    Edge backward(s_, 0, 0.0, next_edge_index);
    e_.push_back(forward);
    e_.push_back(backward);
    outgoing_edges_[s_].push_back(next_edge_index);
    outgoing_edges_[innode_index(ii, 0)].push_back(next_edge_index + 1);
  }

  // add arcs from column c to sink
  for (int ii = 0; ii < r_; ++ii) {
    EdgeIndex next_edge_index = e_.size();
    Edge forward(t_, 1, 0.0, next_edge_index + 1);
    Edge backward(outnode_index(ii, c_ - 1), 0, 0.0, next_edge_index);
    e_.push_back(forward);
    e_.push_back(backward);
    outgoing_edges_[outnode_index(ii, c_ - 1)].push_back(next_edge_index);
    outgoing_edges_[t_].push_back(next_edge_index + 1);
  }

  // add arcs from innodes to outnodes
  node_edges_.resize(r_);
  for (int ii = 0; ii < r_; ++ii) {
    node_edges_[ii].resize(c_);
    for (int jj = 0; jj < c_; ++jj) {
      EdgeIndex next_edge_index = e_.size();
      Edge forward(outnode_index(ii, jj), 1, -abs(a_[ii][jj]),
          next_edge_index + 1);
      Edge backward(innode_index(ii, jj), 0, abs(a_[ii][jj]),
          next_edge_index);
      e_.push_back(forward);
      e_.push_back(backward);
      node_edges_[ii][jj] = next_edge_index;
      outgoing_edges_[innode_index(ii, jj)].push_back(next_edge_index);
      outgoing_edges_[outnode_index(ii, jj)].push_back(next_edge_index + 1);
    }
  }

  // add arcs between columns
  emd_edges_.resize(r_);
  for (int row = 0; row < r_; ++row) {
    emd_edges_[row].resize(c_ - 1);
    for (int col = 0; col < c_ - 1; ++col) {
      size_t ndest = num_destinations(row);
      int first_dest = first_destination(row);
      emd_edges_[row][col].resize(ndest);
      for (size_t idest = 0; idest < ndest; ++idest) {
        EdgeIndex next_edge_index = e_.size();
        Edge forward(innode_index(first_dest + idest, col + 1), 1, 0.0,
            next_edge_index + 1);
        Edge backward(outnode_index(row, col), 0, 0.0, next_edge_index);
        e_.push_back(forward);
        e_.push_back(backward);
        emd_edges_[row][col][idest] = next_edge_index;
        outgoing_edges_[outnode_index(row, col)].push_back(next_edge_index);
        outgoing_edges_[innode_index(first_dest + idest, col + 1)].push_back(
            next_edge_index + 1);
      }
    }
  }

  set_sparsity(0);
}

void EMDFlowNetworkSAP::print_full_graph() {
  printf("Node indices:\n");
  printf("  Source: %lu, sink: %lu\n", s_, t_);
  for (int col = 0; col < c_; ++col) {
    for (int row = 0; row < r_; ++row) {
      printf("  Entry %d,%d: innode: %lu, outnode: %lu\n", row, col,
          innode_index(row, col), outnode_index(row, col));
    }
  }

  printf("Edges:\n");
  for (size_t ii = 0; ii < outgoing_edges_.size(); ++ii) {
    for (size_t jj = 0; jj < outgoing_edges_[ii].size(); ++jj) {
      EdgeIndex curi = outgoing_edges_[ii][jj];
      Edge cur = e_[curi];
      printf("  Edge %lu: from: %lu, to: %lu, cap: %d, cost: %f, "
          "opposite: %lu\n", curi, ii, cur.to, cur.capacity, cur.cost,
          cur.opposite);
    }
  }

  printf("Potentials:\n");
  for (size_t ii = 0; ii < potential_.size(); ++ii) {
    printf("  Node %lu: potential %f\n", ii, potential_[ii]);
  }
}

void EMDFlowNetworkSAP::apply_lambda(double lambda) {
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_ - 1; ++col) {
      size_t ndest = num_destinations(row);
      int first_dest = first_destination(row);
      for (int idest = 0; idest < static_cast<int>(ndest); ++idest) {
        EdgeIndex cur = emd_edges_[row][col][idest];
        e_[cur].cost = lambda * emd_costs_[abs(row - (first_dest + idest))];
        e_[e_[cur].opposite].cost =
            -lambda * emd_costs_[abs(row - (first_dest + idest))];
      }
    }
  }
}

void EMDFlowNetworkSAP::reset_flow() {
  // edges from source to column 1
  for (size_t ii = 0; ii < outgoing_edges_[s_].size(); ++ii) {
    EdgeIndex cur = outgoing_edges_[s_][ii];
    e_[cur].capacity = 1;
    e_[e_[cur].opposite].capacity = 0;
  }

  // edges from column c to sink
  for (size_t ii = 0; ii < outgoing_edges_[t_].size(); ++ii) {
    EdgeIndex cur = outgoing_edges_[t_][ii];
    e_[cur].capacity = 0;
    e_[e_[cur].opposite].capacity = 1;
  }

  // edges from innodes to outnodes
  for (int ii = 0; ii < r_; ++ii) {
    for (int jj = 0; jj < c_; ++jj) {
      EdgeIndex cur = node_edges_[ii][jj];
      e_[cur].capacity = 1;
      e_[e_[cur].opposite].capacity = 0;
    }
  }

  // edges between columns
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_ - 1; ++col) {
      size_t ndest = num_destinations(row);
      for (size_t idest = 0; idest < ndest; ++idest) {
        EdgeIndex cur = emd_edges_[row][col][idest];
        e_[cur].capacity = 1;
        e_[e_[cur].opposite].capacity = 0;
      }
    }
  }
}

void EMDFlowNetworkSAP::compute_initial_potential() {
  // initialize potentials (= distances) to largest possible value
  for (size_t ii = 0; ii < potential_.size(); ++ii) {
    potential_[ii] = numeric_limits<double>::infinity();
  }

  // source and first column have potential 0
  potential_[s_] = 0.0;
  for (int ii = 0; ii < r_; ++ii) {
    potential_[innode_index(ii, 0)] = 0.0;
    potential_[outnode_index(ii, 0)] = -abs(a_[ii][0]);
  }

  // iteratively update next layer based on current layer
  for (int col = 0; col < c_ - 1; ++col) {
    // across column
    for (int row = 0; row < r_; ++row) {
      NodeIndex from = outnode_index(row, col);
      double cur_potential = potential_[from];
      for (size_t ii = 0; ii < outgoing_edges_[from].size(); ++ii) {
        NodeIndex to = e_[outgoing_edges_[from][ii]].to;
        double edge_cost = e_[outgoing_edges_[from][ii]].cost;
        potential_[to] = min(potential_[to], cur_potential + edge_cost);
      }
    }

    // innode to outnode
    for (int row = 0; row < r_; ++row) {
      potential_[outnode_index(row, col + 1)] =
          potential_[innode_index(row, col + 1)] - abs(a_[row][col + 1]);
    }
  }

  // last column to sink
  for (int ii = 0; ii < r_; ++ii) {
    potential_[t_] = min(potential_[t_], potential_[outnode_index(ii, c_ - 1)]);
  }
}

void EMDFlowNetworkSAP::set_sparsity(int s) {
  sparsity_ = s;
}

void EMDFlowNetworkSAP::run_flow(double lambda) {
  typedef pair<double, EMDFlowNetworkSAP::NodeIndex> q_elem;

  reset_flow();
  apply_lambda(lambda);

  //print_full_graph();

  compute_initial_potential();

  vector<EdgeIndex> edge_taken_to(potential_.size(), 0);
  vector<bool> visited(potential_.size(), false);
  vector<double> dst(potential_.size(), numeric_limits<double>::infinity());

  // find a new flow
  for (int total_flow = 0; total_flow < min(sparsity_, r_); ++total_flow) {
    // Dijkstra
    fill(visited.begin(), visited.end(), false);
    fill(dst.begin(), dst.end(), numeric_limits<double>::infinity());
    priority_queue<q_elem> q;

    dst[s_] = 0.0;
    q.push(q_elem(-dst[s_], s_));

    size_t num_found = 0;

    while (!q.empty() && num_found < potential_.size()) {
      q_elem top = q.top();
      q.pop();

      if (visited[top.second]) {
        continue;
      }

      NodeIndex cur_node = top.second;
      visited[cur_node] = true;
      ++num_found;

      NodeIndex next_node;
      for (vector<EdgeIndex>::iterator iter = outgoing_edges_[cur_node].begin();
          iter != outgoing_edges_[cur_node].end(); ++iter) {
        const Edge& e = e_[*iter];
        next_node = e.to;

        ++total_inner_iterations;

        if (e.capacity == 0) {
          continue;
        }
        if (visited[next_node]) {
          continue;
        }
        
        ++checking_inner_iterations;

        double adjusted_edge_cost = e.cost + potential_[cur_node]
            - potential_[next_node];
        if (dst[cur_node] + adjusted_edge_cost < dst[next_node]) {
          dst[next_node] = dst[cur_node] + adjusted_edge_cost;
          q.push(q_elem(-dst[next_node], next_node));
          edge_taken_to[next_node] = *iter;

          ++updating_inner_iterations;
        }
      }
    }

    // change potentials
    for (size_t ii = 0; ii < potential_.size(); ++ii) {
      potential_[ii] += dst[ii];
    }

    // change capacities
    NodeIndex cur_node = t_;
    do {
      Edge& forward_edge = e_[edge_taken_to[cur_node]];
      forward_edge.capacity = 0;
      e_[forward_edge.opposite].capacity = 1;
      cur_node = e_[forward_edge.opposite].to;
    } while (cur_node != s_);
  }

  //print_full_graph();
}

int EMDFlowNetworkSAP::get_EMD_used() {
  int emd_cost = 0;
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_ - 1; ++col) {
      size_t ndest = num_destinations(row);
      int first_dest = first_destination(row);
      for (int idest = 0; idest < static_cast<int>(ndest); ++idest) {
        if (e_[emd_edges_[row][col][idest]].capacity == 0) {
          emd_cost += emd_costs_[abs(row - (first_dest + idest))];
        }
      }
    }
  }
  return emd_cost;
}

double EMDFlowNetworkSAP::get_supported_amplitude_sum() {
  //print_full_graph();

  double amp_sum = 0;
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_; ++col) {
      if (e_[node_edges_[row][col]].capacity == 0) {
        amp_sum += abs(a_[row][col]);
      }
    }
  }
  return amp_sum;
}

void EMDFlowNetworkSAP::get_support(std::vector<std::vector<bool> >* support) {
  if (static_cast<int>(support->size()) != r_) {
    support->resize(r_);
  }
  for (int row = 0; row < r_; ++row) {
    if (static_cast<int>((*support)[row].size()) != c_) {
      (*support)[row].resize(c_);
    }
    for (int col = 0; col < c_; ++col) {
      (*support)[row][col] = (e_[node_edges_[row][col]].capacity == 0);
    }
  }
}

int EMDFlowNetworkSAP::get_num_nodes() {
  return outgoing_edges_.size();
}

int EMDFlowNetworkSAP::get_num_edges() {
  return e_.size();
}

int EMDFlowNetworkSAP::get_num_columns() {
  return c_;
}

int EMDFlowNetworkSAP::get_num_rows() {
  return r_;
}

void EMDFlowNetworkSAP::get_performance_diagnostics(std::string* s) {
  const size_t tmp_size = 2000;
  char tmp[tmp_size];
  snprintf(tmp, tmp_size, "Total inner iterations: %lld\n"
      "Checking inner iterations: %lld\nUpdating inner iterations: %lld\n",
      total_inner_iterations, checking_inner_iterations,
      updating_inner_iterations);
  *s = string(tmp);
}
