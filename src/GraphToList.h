#ifndef GRAPHTOLIST_H
#define GRAPHTOLIST_H

#include <Rcpp.h>
#include <lemon/list_graph.h>
#include <lemon/grid_graph.h>

#include "ProblemData.h"
#include "HelperFunctions.h"

using namespace Rcpp;
using namespace lemon;


/*
 * Exports a lemon graph (for the single-vehicle routing problem) to a list of data frames usable in R.
 */
List routeGraphToList_R(const ProblemData& pd, const std::vector<std::vector<MyGraph::Arc>>& totalPath = std::vector<std::vector<MyGraph::Arc>>());


//
// struct ProblemData{
//     MyGraph& graph_;
//     MyGraph::NodeMap<NodeData>& nodeMap_;
//     MyGraph::ArcMap<ArcData>& arcMap_;
//     MyGraph::Node& startNode_;
//     std::vector<MyGraph::Node> destinations_;
//     int budget_;
//
//     ProblemData(MyGraph& graph,
//                 MyGraph::NodeMap<NodeData>& nodeMap,
//                 MyGraph::ArcMap<ArcData>& arcMap,
//                 MyGraph::Node& startNode,
//                 std::vector<MyGraph::Node> destinations,
//                 int budget)
//         :   graph_(graph), nodeMap_(nodeMap),
//             arcMap_(arcMap),
//             startNode_(startNode),
//             destinations_(destinations), budget_(budget){
//
//     }
//     ProblemData(const ProblemData& pd)
//         :   graph_(pd.graph_), nodeMap_(pd.nodeMap_),
//             arcMap_(pd.arcMap_),
//             startNode_(pd.startNode_), destinations_(pd.destinations_), budget_(pd.budget_){
//
//     }
//
// };



#endif
