/*
 * Given a graph and a permutation of the destinations, calculate the duration of a trip
 */

#include <lemon/grid_graph.h>
#include <vector>
#include <lemon/dijkstra.h>

#include "CalculateLength.h"
#include "GraphToList.h"
#include "HelperFunctions.h"

#include "ProblemData.h"

using namespace lemon;

// the calculated shortest path is to be stored in totalPath
double simulateTrip(const ProblemData& pd,
                    const std::vector<MyGraph::Node>& solution,
                    std::vector<std::vector<MyGraph::Arc>>& totalPath,
                    double& value){

	double result = 0.0;

	totalPath.clear();

	ArcLookUp<MyGraph> arcLookUp(pd.graph_);


	for (std::size_t i = 0; i < solution.size() - 1; i++){
		//std::vector<MyGraph::Node> nodePath;
		std::vector<MyGraph::Arc> currentArcPath;
		double tempResult = calculateShortestPath(pd, solution[i], solution[i+1], currentArcPath);

		if (tempResult != -1){
			result += tempResult;
		} else {
			// return -1;
		}
		// Rcpp::Rcout << "result / tempResult: " << result << " / " << tempResult << std::endl;
		// waitForInput("beginForLoop", DEBUG_ENABLED);
		// for (std::size_t j = 0; j < nodePath.size() - 1; j++){
		//     std::stringstream ss;
		//     ss << "Push back " << j << "/" << nodePath.size() << ": " << pd.graph_.id(nodePath[j]) << " -> " << pd.graph_.id(nodePath[j+1]) << std::endl;
		//     waitForInput(ss.str(), DEBUG_ENABLED);
		//     currentArcPath.push_back(arcLookUp(nodePath[j], nodePath[j+1]));
		//     waitForInput("pushed", DEBUG_ENABLED);
		// }

		// std::vector<MyGraph::Arc> correctedArcPath(currentArcPath.rbegin(), currentArcPath.rend());
		// totalPath.push_back(correctedArcPath);
		totalPath.emplace_back(currentArcPath.rbegin(), currentArcPath.rend());
	}

	value = 0;

	// Maybe skip first and last node of solution (since they correspond to the depot)?
	// for (std::vector<GridGraph::Node>::const_iterator it = solution.begin() + 1; it != solution.end() - 1; it++){
	//     value += nodeMap[*it].value_;
	// }
	for (MyGraph::Node n : solution){
	    value += pd.nodeMap_[n].value_;
	}

	return result;
}

// the calculated shortest path is to be stored in totalPath
double simulateTripDistanceOnly(const ProblemData& pd,
                                const std::vector<MyGraph::Node>& solution,
                                double& value){

	double result = 0.0;

	for (std::size_t i = 0; i < solution.size() - 1; i++){
		double tempResult = calculateShortestPathDistanceOnly(pd, solution[i], solution[i+1]);;
		if (tempResult != -1){
			result += tempResult;
		} else {
			return -1;
		}
	}

	value = 0;
	for (MyGraph::Node n : solution){
	    value += pd.nodeMap_[n].value_;
	}

	return result;
}

/*
 * To-Do: Use a faster algorithm than Dijkstra (e.g., hub labeling)
 */

double calculateShortestPathDistanceOnly(const ProblemData& pd, const MyGraph::Node& source, const MyGraph::Node& target){
    //DijkstraDefaultTraits<MyGraph, MyGraph::ArcMap<double>>::LengthMap lengthMap(pd.graph_);
    DijkstraDouble::LengthMap lengthMap(pd.graph_);

    for (MyGraph::ArcIt ait(pd.graph_); ait != INVALID; ++ait){
        lengthMap[ait] = pd.arcMap_[ait].length_;
    }

    DijkstraDouble d(pd.graph_, lengthMap);
    DijkstraDouble::DistMap dist(pd.graph_);

    d.distMap(dist);
    d.run(source, target);

    return dist[target];
}

double calculateShortestPath(const ProblemData& pd, const MyGraph::Node& source, const MyGraph::Node& target,
                             std::vector<MyGraph::Arc>& resultPath){

    DijkstraDouble::LengthMap lengthMap(pd.graph_);
    //Dijkstra<MyGraph>::LengthMap lengthMap(pd.graph_);
    for (MyGraph::ArcIt ait(pd.graph_); ait != INVALID; ++ait){
        lengthMap[ait] = pd.arcMap_[ait].length_;
        // Rcpp::Rcout << "length: " << lengthMap[ait] << std::endl;
    }

    DijkstraDouble d(pd.graph_, lengthMap);
    DijkstraDouble::DistMap dist(pd.graph_);

    d.distMap(dist);
    d.init();
    d.addSource(source);
    d.start(target);

    // DijkstraDouble::Path dijkstraPath = d.path(target);

    resultPath.clear();
    MyGraph::Arc prevArc = d.predArc(target);
    // waitForInput("path-pushback...", DEBUG_ENABLED);
    // int counter = 0;
    // double actualDistance = 0;

    while (prevArc != INVALID){
        resultPath.push_back(prevArc);
        // Rcpp::Rcout << pd.graph_.id(pd.graph_.source(prevArc)) << " -> " << pd.graph_.id(pd.graph_.target(prevArc)) << std::endl;
        //actualDistance += (int) pd.arcMap_[prevArc].length_;
        prevArc = d.predArc(pd.graph_.source(prevArc));
        // counter++;

        // waitForInput("singleStep", DEBUG_ENABLED);
    }

    // Rcpp::Rcout << "shortest path has " << counter << " arcs and length " << dist[target] << ". "  << std::endl;
    // Rcpp::Rcout << "Path length: " << dijkstraPath.length() << " / " << resultPath.size() << std::endl;
    // waitForInput("path-pushback fertig", DEBUG_ENABLED);

    return dist[target];
}

