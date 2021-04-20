#ifndef CALCULATELENGTH_H
#define CALCULATELENGTH_H

/*
 * Given a graph and a permutation of the destinations, calculate the duration of a trip
 */

#include <lemon/dijkstra.h>
#include <lemon/grid_graph.h>
#include <vector>
#include "GraphToList.h"
#include "ProblemData.h"
#include "HelperFunctions.h"

using namespace lemon;

typedef Dijkstra<MyGraph, MyGraph::ArcMap<double>> DijkstraDouble;

// the calculated shortest path is to be stored in totalPath
double simulateTrip(const ProblemData& pd,
                    const std::vector<MyGraph::Node>& solution,
                    std::vector<std::vector<MyGraph::Arc>>& totalPath,
                    double& value);

double simulateTripDistanceOnly(const ProblemData& pd,
                                const std::vector<MyGraph::Node>& solution,
                                double& value);

double calculateShortestPathDistanceOnly(const ProblemData& pd, const MyGraph::Node& source, const MyGraph::Node& target);

double calculateShortestPath(const ProblemData& pd, const MyGraph::Node& source, const MyGraph::Node& target,
                             std::vector<MyGraph::Arc>& path);

#endif

