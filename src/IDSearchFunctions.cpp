#include "IDSearchFunctions.h"

MyGraph::Node getNodeFromInternalID(const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const int& id){
    for (MyGraph::NodeIt n(graph); n != INVALID; ++n){
        if (nodeMap[n].id_ == id){
            return (MyGraph::Node) n;
        }
    }
    std::stringstream ss;
    ss << "getNodeFrominternalID: Did not find node with id " << id;

    Rcpp::stop(ss.str());
    return graph.nodeFromId(graph.maxNodeId());
}
