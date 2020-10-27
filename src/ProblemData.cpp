#include <Rcpp.h>
#include <list_graph.h>
#include <vector>

#include "ProblemData.h"
#include "HelperFunctions.h"
#include "IO.h"
#include "IDSearchFunctions.h"

using namespace lemon;

void setStartNode(MyGraph::Node& startNode, const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const Rcpp::DataFrame& problemDf){
    Rcpp::IntegerVector valueNodeIDs = problemDf["id"];
    Rcpp::IntegerVector nodeTypeVector = problemDf["type"];
    bool startNodeFound = false;
    int i = 0;

    while(!startNodeFound && i < problemDf.nrow()){
        if (nodeTypeVector(i) == 2) {
            startNode = getNodeFromInternalID(graph, nodeMap, valueNodeIDs(i));//graph.nodeFromId(valueNodeIDs(i));
            startNodeFound = true;
        }
        i++;
    }
    // for (MyGraph::NodeIt n(graph); n != INVALID; ++n){
    //     Rcpp::Rcout << graph.id(n) << " ";
    // }
    if (!startNodeFound){
        Rcpp::stop("start node was not found");
    }
}

// destination do not contain the starting node
/*
 * Why? I dunno anymore...
 */
void fillDestinations(std::vector<MyGraph::Node>& destinations, const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const Rcpp::DataFrame& problemDf){
    Rcpp::IntegerVector valueNodeIDs = problemDf["id"];
    for (int i = 0; i < valueNodeIDs.size(); i++){
        MyGraph::Node n = getNodeFromInternalID(graph, nodeMap, valueNodeIDs(i));
        if (nodeMap[n].type_ == 1){
            destinations.push_back(getNodeFromInternalID(graph, nodeMap, valueNodeIDs(i)));
        }
    }
}

void readDynamicChangesFromFile(ProblemData& problemData, std::string path){

    if (path != ""){
        std::string line;
        std::ifstream infile(path);
        std::getline(infile, line);
        std::vector<int> counterForChanges = {0,0,0,0};
        while (std::getline(infile, line)) {
            stringToChangeConverter conv;
            // Rcpp::Rcout << line << std::endl;
            Change change = conv(line);
            // Rcpp::Rcout << change << std::endl;

            /*
             * if a change is set for multiple time measures, the (arbitrarily set) preference is:
             * time -> evaluation -> bit vector
             */
            if (change.atTime_ != 0){
                problemData.dynamicChangesByTime_.push_back(change);
                counterForChanges[0]++;
            } else {
                if (change.atEvaluation_ != 0){
                    problemData.dynamicChangesByEvaluation_.push_back(change);
                    counterForChanges[1]++;
                } else {
                    if (change.atTestedBitVector_ != 0) {
                        problemData.dynamicChangesByTestedBitVectors_.push_back(change);
                        counterForChanges[2]++;
                    } else {
                        if (change.atDistanceEvaluation_ != 0) {
                            problemData.dynamicChangesByDistanceEvaluation_.push_back(change);
                            counterForChanges[3]++;
                        } else {
                            Rcpp::Rcerr << "WARNING: A change in the change file has no valid starting time and is thus ignored:" << std::endl;
                            Rcpp::Rcerr << line << std::endl;
                        }
                    }
                }
            }
            // bool changeIsMixed = true;
            //
            // if (change.atEvaluation_ == 0 || change.atTestedBitVector_ == 0){
            //     problemData.dynamicChangesByTime_.push_back(change);
            //     changeIsMixed = false;
            // }
            // if (change.atTime_ == 0 || change.atTestedBitVector_ == 0){
            //     problemData.dynamicChangesByEvaluation_.push_back(change);
            //     changeIsMixed = false;
            // }
            //
            // if (changeIsMixed){
            //     problemData.dynamicChangesMixed_.push_back(change);
            // }
        }

        // make termination criterion dependent on the time measure
        // that has the most changes
        int maxNumber = 0;
        int maxIndex = 0;
        for (std::size_t d = 0; d < counterForChanges.size(); d++){
            if (counterForChanges[d] > maxNumber) {
                maxNumber = counterForChanges[d];
                maxIndex = d;
            }
        }
        switch(maxIndex) {
        case(1):
            problemData.terminationCriterion_ = "evaluation";
            break;
        case(2):
            problemData.terminationCriterion_ = "bitVector";
            break;
        case(3):
            problemData.terminationCriterion_ = "distanceEvaluation";
            break;
        default:
            problemData.terminationCriterion_ = "combined";
        }
        Rcpp::Rcout << "terminationCriterion: " << problemData.terminationCriterion_ << std::endl;

    }
}
