// #include <Rcpp.h>
// #include <grid_graph.h>
// #include <list_graph.h>
// #include <stdlib.h>
// #include <lgf_writer.h>
// #include <lgf_reader.h>
// #include <algorithm>
// #include <random>
//
// #include <chrono>
// #include <thread>
//
// #include "GraphToList.h"
// #include "IO.h"
// #include "HelperFunctions.h"
//
// #include "CalculateLength.h"
// #include "ShortestPath.h"
// #include "SimulatePath.h"
//
// using namespace Rcpp;
// using namespace lemon;
//
//
//
// void simulatePath(ListDigraph& graph,
//                   ListDigraph::NodeMap<NodeData>& nodeMap,
//                   ListDigraph::ArcMap<ArcData>& arcMap,
//                   std::vector<ListDigraph::Node>& destinations){
//
//     long numberOfEvaluation = 0;
//     std::vector<std::vector<ListDigraph::Arc>> totalPath;
//     double value = 0;
//
//     int stepCounter = 0;
//     int timeUnits = 0;
//     bool destinationReached = false;
//     // long tickLength = 10;
//     double totalValue = 0;
//
//     double tempLength = 0;
//
//     Rcout << "Plot..." << std::endl;
//
//     ProblemData pd(graph, nodeMap, arcMap, destinations[0], destinations, -1);
//
//
//     tempLength = simulateTrip(pd, destinations, totalPath, value);
//     for (std::size_t i = 0; i < totalPath.size(); i++){
//         // for (ListDigraph::Arc a : totalPath[i]){
//         //     arcMap[a].isOnRouteOf_ = i+1;
//         // }
//     }
//
//     // plotGridRoute_cpp(graph, nodeMap, arcMap, forbiddenEdges);
//     // waitForInput("startSimulation", DEBUG_ENABLED);
// //
// //     while (!destinationReached){
// //
// //         // mark all destinations (perhaps not necessary)
// //         for (std::size_t i = 0; i < destinations.size(); i++){
// //             nodeMap[destinations[i]].type_ = 1;
// //         }
// //
// //         destinations[0] = graph.target(totalPath[0][0]);
// //         nodeMap[destinations[0]].type_ = 2;
// //         arcMap[totalPath[0][0]].isOnRouteOf_ = 0;
// //
// //         stepCounter++;
// //
// //         if (destinations[0] == destinations[1]){
// //             // final destination
// //             if (destinations[0] == destinations.back()){
// //                 destinationReached = true;
// //             }
// //             //nodeMap[destinations[1]].type_ = 0;
// //             totalValue += nodeMap[destinations[1]].value_;
// //             nodeMap[destinations[1]].value_ = 0;
// //             destinations.erase(destinations.begin() + 1);
// //         }
// //
// //
// //         tempLength = simulateTrip(graph, nodeMap, forbiddenEdges, destinations, totalPath, value);
// //
// //         if (tempLength != -1){
// //             // construct path from sequence of arcs
// //             for (std::size_t i = 0; i < totalPath.size(); i++){
// //                 for (GridGraph::Arc a : totalPath[i]){
// //                     arcMap[a].isOnRouteOf_ = i+1;
// //                 }
// //             }
// //
// //             // plot
// //             plotGridRoute_cpp(graph, nodeMap, arcMap, forbiddenEdges);
// //             std::this_thread::sleep_for(std::chrono::milliseconds(200));
// //         }
// //
// //         timeUnits++;
// //         numberOfEvaluation++;
// //     }
// //     Rcout << "evals - ticks - value: " << numberOfEvaluation << " - " << timeUnits  << " - " << totalValue << std::endl;
//
// }
