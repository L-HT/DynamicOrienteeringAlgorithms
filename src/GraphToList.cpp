/*
 * Functions to convert a lemon graph to an R list.
 * This list contains 3 elements:
 * 	-GraphData (metadata about the graph)
 * 	-divergenceDf (which nodes a sources/sinks, which nodes are true supply nodes)
 * 	-arcDf (capacity and current flow for each arc
 *
 * This data is later needed to plot the network
 */

#include <Rcpp.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <grid_graph.h>

#include "GraphToList.h"
#include "HelperFunctions.h"
#include "ProblemData.h"
#include "IO.h"

using namespace Rcpp;
using namespace lemon;

List routeGraphToList_R(const ProblemData& pd, const std::vector<std::vector<MyGraph::Arc>>& totalPath){

	int numberOfValidNodes = 0;

	for (MyGraph::NodeIt n(pd.graph_); n != INVALID; ++n){
	    numberOfValidNodes++;
    }

	// node data
	std::vector<int> nodeIds = std::vector<int>(numberOfValidNodes,0);
	std::vector<double> xPositions = std::vector<double>(numberOfValidNodes,0.0);
	std::vector<double> yPositions = std::vector<double>(numberOfValidNodes,0.0);
	std::vector<int> typeVector = std::vector<int>(numberOfValidNodes, 0);
	std::vector<int> valueVector = std::vector<int>(numberOfValidNodes, 0);
	std::size_t counter = 0;

	// waitForInput("startNodeData", true);

	for (MyGraph::NodeIt n(pd.graph_); n != INVALID; ++n){
		nodeIds[counter] = pd.nodeMap_[n].id_;//pd.graph_.id(n);
		NodeData nd = pd.nodeMap_[n];

		yPositions[counter] = nd.coordinates_.y;
		xPositions[counter] = nd.coordinates_.x;

		typeVector[counter] = nd.type_;
		valueVector[counter] = nd.value_;
		counter++;

	}


	DataFrame nodeDf = DataFrame::create(
		Named("nodeId") = nodeIds,
		Named("x") = xPositions,
		Named("y") = yPositions,
		Named("type") = typeVector,
		Named("value") = valueVector
	);

	/////////////// arc data

	int numberOfArcs = 0;

	for (MyGraph::ArcIt a(pd.graph_); a != INVALID; ++a){
	    numberOfArcs++;
	}

	std::vector<int> starts = std::vector<int>(numberOfArcs, 0);
	std::vector<int> ends = std::vector<int>(numberOfArcs, 0);
	std::vector<double> lengths = std::vector<double>(numberOfArcs, 0.0);
	std::vector<int> onRouteOfVector = std::vector<int>(numberOfArcs, 0);

	counter = 0;
	// waitForInput("fill arcmap...", DEBUG_ENABLED);
	MyGraph::ArcMap<int> onRouteOfArcMap(pd.graph_, 0);
	if (totalPath.size() != 0){
    	for (auto outer = totalPath.begin(); outer != totalPath.end(); outer++){
    	    for (auto inner = (*outer).begin(); inner != (*outer).end(); inner++){
    	        // Rcpp::Rcout << counter << ": " << pd.graph_.id(pd.graph_.source(*inner)) << " -> " << pd.graph_.id(pd.graph_.target(*inner)) << std::endl;
    	        onRouteOfArcMap[*inner] = ++counter;
    	    }
    	}
	}
	// waitForInput("filled arcmap...", DEBUG_ENABLED);

	counter = 0;
	for (MyGraph::ArcIt a(pd.graph_); a != INVALID; ++a){

	//	ArcData ad = pd.arcMap_[a];
		starts[counter] = pd.arcMap_[a].fromID_;
		ends[counter] = pd.arcMap_[a].toID_;
		lengths[counter] = pd.arcMap_[a].length_;
        onRouteOfVector[counter] = onRouteOfArcMap[a];

		counter++;
	}

	// resize
	// arcIds.resize(counter);
	// starts.resize(counter);
	// ends.resize(counter);
	//
	// row1.resize(counter);
	// col1.resize(counter);
	// row2.resize(counter);
	// col2.resize(counter);

	//existVector.resize(counter);
	onRouteOfVector.resize(counter);

	DataFrame arcDf = DataFrame::create(
		Named("from") = starts,
		Named("to") = ends,
		Named("length") = lengths,
		Named("onRouteOf") = onRouteOfVector
	);

	// waitForInput("finishArcData", true);

	return List::create(
		Named("nodeDf") = nodeDf,
		Named("arcDf") = arcDf
	);

}

// void plotNetworkRoute_cpp(const ProblemData& pd, const std::vector<std::vector<MyGraph::Arc>>& totalPath){
// 	Environment e(MYPACKAGE);
// 	Function plotNetworkRoute = e[PLOTFUNCTION];
// 	// waitForInput("list1", DEBUG_ENABLED);
// 	List graphAsList = routeGraphToList_R(pd, totalPath);
// 	// waitForInput("list2", DEBUG_ENABLED);
// 	DataFrame nodeDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["nodeDf"]);
// 	DataFrame arcDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["arcDf"]);
// 	plotNetworkRoute(R_NilValue, nodeDf, arcDf);
// }
//
//
//
// void plotNetworkRoute_R(const Rcpp::DataFrame& nodeDf,
//                         const Rcpp::DataFrame& arcDf,
//                         const Rcpp::DataFrame& problemDf,
//                         const std::vector<int>& nodeIDs){
//     MyGraph graph;
//     MyGraph::ArcMap<ArcData> arcMap(graph);
//     MyGraph::NodeMap<NodeData> nodeMap(graph);
//     MyGraph::Node startNode;
//     std::vector<MyGraph::Node> destinations;
//     readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);
//     readProblemDataIntoMaps(problemDf, graph, nodeMap, arcMap);
//     fillDestinations(destinations, graph, nodeMap, problemDf);
//     setStartNode(startNode, graph, nodeMap, problemDf);
//
//     ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, 0, "");
//
//     Environment e(MYPACKAGE);
//     Function plotNetworkRoute = e[PLOTFUNCTION];
//
//     std::vector<MyGraph::Node> mySolution;
//     //mySolution.push_back(pd.startNode_);
//     for (int i : nodeIDs){
//         mySolution.push_back(pd.graph_.nodeFromId(i));
//     }
//     //mySolution.push_back(pd.startNode_);
//
//     std::vector<std::vector<MyGraph::Arc>> totalPath;
//     double value = 0;
//     double length = simulateTrip(pd, mySolution, totalPath, value);
//     Rcpp::Rcout << "length/value: " << length << "/" << value << std::endl;
//     writeSolution(pd, mySolution, "./solutions/ii");
//     plotNetworkRoute_cpp(pd, totalPath);
//
//     // waitForInput("list1", DEBUG_ENABLED);
//     List graphAsList = routeGraphToList_R(pd, totalPath);
//     // waitForInput("list2", DEBUG_ENABLED);
//     DataFrame nodeDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["nodeDf"]);
//     DataFrame arcDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["arcDf"]);
//     plotNetworkRoute(R_NilValue, nodeDf, arcDf);
// }
//
