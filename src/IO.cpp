#include <Rcpp.h>
#include <grid_graph.h>
#include <list_graph.h>
#include <stdlib.h>
#include <lgf_writer.h>
#include <lgf_reader.h>

#include "GraphToList.h"
#include "HelperFunctions.h"

#include "IO.h"
#include "IDSearchFunctions.h"

using namespace Rcpp;
using namespace lemon;


ProblemData problemDataFromLGF(std::string filePath, MyGraph& graph, MyGraph::ArcMap<ArcData>& arcMap,
                               MyGraph::NodeMap<NodeData>& nodeMap, MyGraph::Node& startNode,
                               std::vector<MyGraph::Node>& destinations){



    double budget = 0;
    // waitForInput("call digraphreader", true);
    digraphReader(graph, filePath)
        .nodeMap("nodeData", nodeMap, stringToNodeDataConverter())
        .arcMap("arcData", arcMap, stringToArcDataConverter())
        .attribute("budget", budget)
        // .attribute("gridCols", cols)
        .run();


    // waitForInput("search dests and start", true);
    int counter = 0;
    for (MyGraph::NodeIt it(graph); it != INVALID; ++it){
        if (nodeMap[it].type_ == 2){
            startNode = it;
        }
        if (nodeMap[it].type_ == 1){
            destinations.push_back(it);
        }
        counter++;
    }

    Rcpp::Rcout << "graph hat " << counter << " Knoten" << std::endl;

    // waitForInput("get problem name", true);
    std::string problemName = filePath;
    int lastSlashPos = problemName.find_last_of("\\/");
    if (lastSlashPos != std::string::npos) {
        problemName.erase(0, lastSlashPos + 1);
    }

    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, budget, problemName);

    // waitForInput("und...", true);
    // for (MyGraph::NodeIt it(pd.graph_); it != INVALID; ++it){
    //     counter++;
    // }
    // Rcpp::Rcout << "graph hat " << counter << " Knoten" << std::endl;
    // waitForInput("tschuessss.", true);
    return pd;
}

//' @export
// [[Rcpp::export]]
List convertLGFtoR(std::string filePath){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    ProblemData pd = problemDataFromLGF(filePath, graph, arcMap, nodeMap, startNode, destinations);

    // int counter = 0;
    // for (MyGraph::NodeIt it(pd.graph_); it != INVALID; ++it){
    //     counter++;
    // }
    // Rcpp::Rcout << "graph hat " << counter << " Knoten" << std::endl;

    // waitForInput("pd gelesen", true);
    std::vector<int> nodeIDs, problemIDs, fromIDs, toIDs, typeValues;
    std::vector<double> xValues, yValues, nodeValues, lengthValues;
    for (MyGraph::NodeIt it(pd.graph_); it != INVALID; ++it){
        // Rcpp::Rcout << pd.nodeMap_[it] << std::endl;
        // waitForInput("", true);
        if (pd.nodeMap_[it].type_ >= 1){
            problemIDs.push_back(pd.nodeMap_[it].id_);
            nodeValues.push_back(pd.nodeMap_[it].value_);
            typeValues.push_back(pd.nodeMap_[it].type_);
        }
        nodeIDs.push_back(pd.nodeMap_[it].id_);
        xValues.push_back(pd.nodeMap_[it].coordinates_.x);
        yValues.push_back(pd.nodeMap_[it].coordinates_.y);
    }

    // waitForInput("zweite for-Schleife", true);
    for (MyGraph::ArcIt it(pd.graph_); it != INVALID; ++it){
        fromIDs.push_back(pd.arcMap_[it].fromID_);
        toIDs.push_back(pd.arcMap_[it].toID_);
        lengthValues.push_back(pd.arcMap_[it].length_);
    }

    // waitForInput("arc df", true);
    DataFrame arcDf = DataFrame::create(
        Named("from") = fromIDs,
        Named("to") = toIDs,
        Named("length") = lengthValues
    );

    // waitForInput("node df", true);
    DataFrame nodeDf = DataFrame::create(
        Named("id") = nodeIDs,
        Named("x") = xValues,
        Named("y") = yValues
    );

    // waitForInput("problemDf", true);
    DataFrame problemDf = DataFrame::create(
        Named("id") = problemIDs,
        Named("type") = typeValues,
        Named("value") = nodeValues
    );
    return List::create(
        Named("nodeDf") = nodeDf,
        Named("arcDf") = arcDf,
        Named("problemDf") = problemDf,
        Named("budget") = pd.budget_
    );
}

/*
 * reads an instance into c++ data structures
 */
void readGridGraphFile(std::string filePath, GridGraph& graph,
                       GridGraph::NodeMap<NodeRouteData>& nodeMap,
                       GridGraph::ArcMap<ArcRouteData>& arcMap,
                       GridGraph::NodeMap<bool>& forbiddenNodes,
                       GridGraph::ArcMap<int>& forbiddenEdges,
                       GridGraph::Node& startNode)
                       {
	int rows = 10;
	int cols = 10;
	MyGraph graph_temp;

	MyGraph::NodeMap<int> rowCoordinates_temp(graph_temp);
	MyGraph::NodeMap<int> colCoordinates_temp(graph_temp);
	MyGraph::NodeMap<int> nodeType_temp(graph_temp);
	MyGraph::NodeMap<int> nodeValue_temp(graph_temp);

	MyGraph::NodeMap<NodeRouteData> nodeMap_temp(graph_temp);
	MyGraph::ArcMap<ArcRouteData> arcMap_temp(graph_temp);

	MyGraph::ArcMap<int> forbiddenEdges_temp(graph_temp, 0);
	MyGraph::NodeMap<bool> forbiddenNodes_temp(graph_temp, false);

	MyGraph::ArcMap<int> arcRow1_temp(graph_temp);
	MyGraph::ArcMap<int> arcCol1_temp(graph_temp);
	MyGraph::ArcMap<int> arcRow2_temp(graph_temp);
	MyGraph::ArcMap<int> arcCol2_temp(graph_temp);

	digraphReader(graph_temp, filePath)
		  .nodeMap("row", rowCoordinates_temp)
		  .nodeMap("col", colCoordinates_temp)
		  .nodeMap("type", nodeType_temp)
          .nodeMap("value", nodeValue_temp)
		  .nodeMap("forbidden", forbiddenNodes_temp)
		  .arcMap("row1", arcRow1_temp)
		  .arcMap("col1", arcCol1_temp)
		  .arcMap("row2", arcRow2_temp)
		  .arcMap("col2", arcCol2_temp)
		  .arcMap("forbidden", forbiddenEdges_temp)
		  .attribute("gridRows", rows)
		  .attribute("gridCols", cols)
		  .run();

	//////////////////////////////////

	// due to the nature of lemon::GridGraph, it is necessary to copy
	// everything from graph_temp by hand

	graph.resize(rows, cols);

	GridGraph::Node sourceNode;
	std::vector<GridGraph::Node> targetNodes;

	for (MyGraph::NodeIt it(graph_temp); it != INVALID; ++it){
		GridGraph::Node correspondingNode = graph(colCoordinates_temp[it], rowCoordinates_temp[it]);
		nodeMap[correspondingNode].coordinates_.y = rowCoordinates_temp[it];
		nodeMap[correspondingNode].coordinates_.x = colCoordinates_temp[it];
		nodeMap[correspondingNode].type_ = nodeType_temp[it];
		nodeMap[correspondingNode].value_ = nodeValue_temp[it];

		forbiddenNodes[correspondingNode] = forbiddenNodes_temp[it];

		// start node
		if (nodeType_temp[it] == 2){
			startNode = correspondingNode;
		}
	}

	// determining the corresponding arc just by id seems to be valid
	// (even though this has not been thoroughly checked so far)
	for (MyGraph::ArcIt it(graph_temp); it != INVALID; ++it){

		GridGraph::Arc correspondingArc = graph.arcFromId(graph_temp.id(it));
		arcMap[correspondingArc].col1_ = arcCol1_temp[it];
		arcMap[correspondingArc].col2_ = arcCol2_temp[it];

		arcMap[correspondingArc].row1_ = arcRow1_temp[it];
		arcMap[correspondingArc].row2_ = arcRow2_temp[it];

		forbiddenEdges[correspondingArc] = forbiddenEdges_temp[it];
	}
}



//' @export
// [[Rcpp::export]]
void saveDigraphAsLGF(std::string filePath, Rcpp::DataFrame& nodeDf, Rcpp::DataFrame& arcDf){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);

    digraphWriter(graph, filePath)
        .nodeMap("nodeData", nodeMap, nodeArcDataToStringConverter())
        .arcMap("arcData", arcMap, nodeArcDataToStringConverter())
        .run();
}

//' @export
// [[Rcpp::export]]
void saveDigraphWithProblemDataAsLGF(std::string filePath, Rcpp::DataFrame& nodeDf, Rcpp::DataFrame& arcDf, Rcpp::DataFrame& problemDf, double budget){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);
    readProblemDataIntoMaps(problemDf,
                            graph,
                            nodeMap,
                            arcMap);

    digraphWriter(graph, filePath)
        .nodeMap("nodeData", nodeMap, nodeArcDataToStringConverter())
        .arcMap("arcData", arcMap, nodeArcDataToStringConverter())
        .attribute("budget", budget)
        .run();
}


/*
 * For a given nodeDf and arcDf (from R), construct a MyGraph containing this data.
 */
void readMyGraphIntoVariables(const Rcpp::DataFrame& nodeDf,
                                  const Rcpp::DataFrame& arcDf,
                                  MyGraph& graph,
                                  MyGraph::NodeMap<NodeData>& nodeMap,
                                  MyGraph::ArcMap<ArcData>& arcMap){

    // MyGraph graph;
    // MyGraph::ArcMap<ArcData> arcMap(graph);
    // MyGraph::NodeMap<NodeData> nodeMap(graph);

    Rcpp::IntegerVector ids = nodeDf["id"];
    Rcpp::NumericVector xValues = nodeDf["x"];
    Rcpp::NumericVector yValues = nodeDf["y"];
    // waitForInput("now nodeDf", DEBUG_ENABLED);

    for (int i = 0; i < nodeDf.nrow(); i++){
        MyGraph::Node n = graph.addNode();

        nodeMap[n] = NodeData(ids(i), xValues(i), yValues(i), 0, 0);//NodeData(graph.id(n), xValues(i), yValues(i), 0, 0);

    }

    // waitForInput("now arcDf", DEBUG_ENABLED);
    Rcpp::IntegerVector fromIDs = arcDf["from"];
    Rcpp::IntegerVector toIDs = arcDf["to"];
    Rcpp::NumericVector distanceValues = arcDf["length"];

    for (int i = 0; i < arcDf.nrow(); i++){
        // Rcpp::Rcout << "from - to - length: " << fromIDs(i) << " - " << toIDs(i) << " - " << distanceValues(i) << std::endl;
        MyGraph::Node n1 = getNodeFromInternalID(graph, nodeMap, fromIDs(i));//graph.nodeFromId(fromIDs(i));//getNodeFromInternalID(graph, nodeMap, fromIDs(i));
        MyGraph::Node n2 = getNodeFromInternalID(graph, nodeMap, toIDs(i));//graph.nodeFromId(toIDs(i));//getNodeFromInternalID(graph, nodeMap, toIDs(i));
        if (n1 == INVALID){
            Rcpp::Rcerr << "node n1 with fromID " << fromIDs(i) << " is invalid.";
            Rcpp::stop("n1 invalid");
        }
        if (n2 == INVALID){
            Rcpp::Rcerr << "node n2 with fromID " << fromIDs(i) << " is invalid.";
            Rcpp::stop("n2 invalid");
        }
        MyGraph::Arc a = graph.addArc(n1, n2);

        arcMap[a] = ArcData(fromIDs(i), toIDs(i), distanceValues(i));
    }

}

// for a given distance matrix, construct a graph that only contains the
// nodes from the distance matrix
void readDistanceMatrixIntoVariables(const Rcpp::DataFrame& nodeDf,
                              const Rcpp::DataFrame& arcDf,
                              const Rcpp::DataFrame& problemDf,
                              MyGraph& graph,
                              MyGraph::NodeMap<NodeData>& nodeMap,
                              MyGraph::ArcMap<ArcData>& arcMap,
                              const Rcpp::NumericMatrix& distanceMatrix){


    Rcpp::IntegerVector ids = nodeDf["id"];
    Rcpp::NumericVector xValues = nodeDf["x"];
    Rcpp::NumericVector yValues = nodeDf["y"];

    Rcpp::IntegerVector problemDf_ids = problemDf["id"];
    Rcpp::NumericVector typeValues = problemDf["type"];
    Rcpp::NumericVector valueValues = problemDf["value"];

    for (int i = 0; i < distanceMatrix.nrow(); i++){
        MyGraph::Node n = graph.addNode();

        // search problemDf-node in nodeDf
        Rcpp::IntegerVector::const_iterator it = std::find(ids.begin(), ids.end(), problemDf_ids(i));
        if (it == ids.end()) {
            std::stringstream ss;
            ss << "readDistanceMatrixIntoVariables: nodeDf does not contain node with ID " << problemDf_ids(i) << "!\n";
            Rcpp::stop(ss.str());
        }

        std::size_t j = it - ids.begin();
        nodeMap[n] = NodeData(problemDf_ids(i), xValues(j), yValues(j), valueValues(i), typeValues(i));
        // Rcpp::Rcout << "added: " << nodeMap[n] << "\n";
    }

    for (int i = 0; i < distanceMatrix.nrow(); i++) {
        for (int j = 0; j < distanceMatrix.ncol(); j++) {
            // Rcpp::Rcout << "from - to - length: " << fromIDs(i) << " - " << toIDs(i) << " - " << distanceValues(i) << std::endl;
            MyGraph::Node n1 = getNodeFromInternalID(graph, nodeMap, problemDf_ids(i));//graph.nodeFromId(fromIDs(i));//getNodeFromInternalID(graph, nodeMap, fromIDs(i));
            MyGraph::Node n2 = getNodeFromInternalID(graph, nodeMap, problemDf_ids(j));//graph.nodeFromId(toIDs(i));//getNodeFromInternalID(graph, nodeMap, toIDs(i));
            if (n1 == INVALID){
                Rcpp::Rcerr << "node n1 with fromID " << problemDf_ids(i) << " is invalid.";
                Rcpp::stop("n1 invalid");
            }
            if (n2 == INVALID){
                Rcpp::Rcerr << "node n2 with fromID " << problemDf_ids(j) << " is invalid.";
                Rcpp::stop("n2 invalid");
            }
            MyGraph::Arc a = graph.addArc(n1, n2);

            arcMap[a] = ArcData(problemDf_ids(i), problemDf_ids(j), distanceMatrix(i,j));
        }
    }
}




/*
* For a given nodeDf and arcDf (from R), construct a MyGraph containing this data.
*/
void readProblemDataIntoMaps(const Rcpp::DataFrame& problemDf,
                                  MyGraph& graph,
                                  MyGraph::NodeMap<NodeData>& nodeMap,
                                  MyGraph::ArcMap<ArcData>& arcMap){

    Rcpp::IntegerVector ids = problemDf["id"];
    Rcpp::NumericVector typeValues = problemDf["type"];
    Rcpp::NumericVector valueValues = problemDf["value"];

    for (int i = 0; i < problemDf.nrow(); i++){
        MyGraph::Node n = getNodeFromInternalID(graph, nodeMap, ids(i));//graph.nodeFromId(ids(i));
        nodeMap[n].id_ = ids(i);
        nodeMap[n].type_ = typeValues(i);
        nodeMap[n].value_ = valueValues(i);
    }

    /*
     * arcMap is not modified so far...
     */

}
