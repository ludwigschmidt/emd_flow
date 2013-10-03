#include "emd_flow_network_factory.h"
#include "emd_flow_network.h"
#include "emd_flow_network_sap.h"

#include <memory>

#ifdef USE_LEMON
  #include "emd_flow_network_lemon.h"
  #include <lemon/network_simplex.h>
  #include <lemon/cost_scaling.h>
  #include <lemon/capacity_scaling.h>
  using namespace lemon;
#endif

using namespace std;

auto_ptr<EMDFlowNetwork> EMDFlowNetworkFactory::create_EMD_flow_network(
        const vector<vector<double> >& amplitudes,
        int outdegree_vertical_distance,
        const std::vector<double>& emd_costs,
        EMDFlowNetworkType type) {
  #ifdef USE_LEMON
    if (type == kLemonCostScaling) {
      return auto_ptr<EMDFlowNetwork>(
          new EMDFlowNetworkLemon<CostScaling<ListDigraph, int, double> >(
              amplitudes, outdegree_vertical_distance, emd_costs));
    } else if (type == kLemonNetworkSimplex) {
      return auto_ptr<EMDFlowNetwork>(
          new EMDFlowNetworkLemon<NetworkSimplex<ListDigraph, int, double> >(
              amplitudes, outdegree_vertical_distance, emd_costs));
    } else if (type == kLemonCapacityScaling) {
      return auto_ptr<EMDFlowNetwork>(
          new EMDFlowNetworkLemon<CapacityScaling<ListDigraph, int, double> >(
              amplitudes, outdegree_vertical_distance, emd_costs));
    } else
  #endif

  if (type == kShortestAugmentingPath) {
    return auto_ptr<EMDFlowNetwork>(new EMDFlowNetworkSAP(amplitudes,
        outdegree_vertical_distance, emd_costs));
  } else {
    return auto_ptr<EMDFlowNetwork>();
  }
}

EMDFlowNetworkFactory::EMDFlowNetworkType EMDFlowNetworkFactory::parse_type(
    const std::string& name) {
  if (name == "lemon-costscaling") {
    return kLemonCostScaling;
  } else if (name == "lemon-networksimplex") {
    return kLemonNetworkSimplex;
  } else if (name == "lemon-capacityscaling") {
    return kLemonCapacityScaling;
  } else if (name == "sap" || name == "shortest-augmenting-path") {
    return kShortestAugmentingPath;
  } else {
    return kUnknownType;
  }
}

