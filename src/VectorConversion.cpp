#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "ProblemData.h"

// convert solution to vector with indices referring to destinations vector (containing the starting node as first element)
std::vector<int> nodeVectorToIndexVectorWithLimit(const ProblemData& pd, const std::vector<MyGraph::Node>& solution, std::size_t upperLimit){
    std::vector<int> result(solution.size(), 0);
    int outerCounter = 0;

    //for (MyGraph::Node n : solution){
    for (std::size_t i = 0; i < upperLimit; i++){
        MyGraph::Node n = solution[i]; // maybe a pointer here?
        bool nodeFound = false;
        int innerCounter = 0;
        do{
            // Rcpp::Rcout << "innerCounter: " << innerCounter << " / " << pd.destinations_.size() << std::endl;
            if (n == pd.destinations_[innerCounter]){
                // + 1 because index 0 is used for the starting node (which is not contained in pd.destinations_))
                result[outerCounter++] = innerCounter + 1;
                // Rcpp::Rcout << "write: " << innerCounter + 1 << std::endl;
                nodeFound = true;
            } else {
                innerCounter++;
                if (innerCounter >= (int) pd.destinations_.size()){
                    // if (pd.graph_.id(n) == pd.graph_.id(pd.startNode_)){
                    if (pd.nodeMap_[n].id_ == pd.nodeMap_[pd.startNode_].id_){
                        result[outerCounter++] = 0;
                        // Rcpp::Rcout << "write: " << 0 << std::endl;
                        nodeFound = true;
                    } else {
                        std::stringstream ss;
                        ss << "Error in VectorConversion.cpp: Cannot find node ID " << pd.nodeMap_[n].id_ << " in destionations vector! (it's also not the start node)" << std::endl;
                        Rcpp::Rcerr << "solution: ";
                        for (MyGraph::Node n2 : solution){
                            Rcpp::Rcerr << pd.nodeMap_[n2].id_ << " ";
                        }
                        Rcpp::Rcerr << "\ndestinations:\n";
                        for (MyGraph::Node n2 : pd.destinations_){
                            Rcpp::Rcerr << pd.nodeMap_[n2].id_ << " ";
                        }
                        Rcpp::Rcerr << "\nstart node:\n" << pd.nodeMap_[pd.startNode_].id_;
                        Rcpp::Rcerr << std::endl;
                        Rcpp::stop(ss.str());
                    }
                }
            }
        } while (!nodeFound);
    }

    return result;
}

std::vector<int> nodeVectorToIndexVectorWithStartNode(const ProblemData& pd, const std::vector<MyGraph::Node>& solution){
    return nodeVectorToIndexVectorWithLimit(pd, solution, solution.size() - 1);
}
std::vector<int> nodeVectorToIndexVector(const ProblemData& pd, const std::vector<MyGraph::Node>& solution){
    return nodeVectorToIndexVectorWithLimit(pd, solution, solution.size());
}

std::vector<MyGraph::Node> indexVectorToNodeVector(const ProblemData& pd, const int* indices,
                                                   const std::size_t numberOfElements,
                                                   const std::map<int, int>& indicesMap){

    std::vector<MyGraph::Node> result;

    // shift the solution such that the depot ist the first node in the solution
    std::size_t startPos = 0;
    for (std::size_t i = 0; i < numberOfElements; i++){
        if (indices[i] == 0){
            startPos = i;
        }
    }

    // Rcpp::Rcout << "startPos: " << startPos << ", " << numberOfElements << " elements." << std::endl;

    result.push_back(pd.startNode_);
    for (std::size_t i = startPos+1; i < numberOfElements; i++){
        int correspondingIndexInDestinations = indicesMap.at(indices[i]) - 1;

        result.push_back(pd.destinations_[correspondingIndexInDestinations]);
        // Rcpp::Rcout << pd.graph_.id(pd.destinations_[correspondingIndexInDestinations]) << " - " << pd.graph_.id(result.back()) << std::endl;
    }
    for (std::size_t i = 0; i < startPos; i++){
        int correspondingIndexInDestinations = indicesMap.at(indices[i]) - 1;
        result.push_back(pd.destinations_[correspondingIndexInDestinations]);
    }
    result.push_back(pd.startNode_);
    return result;
}

