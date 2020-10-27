#ifndef VECTOR_CONVERSION_H
#define VECTOR_CONVERSION_H

#include "ProblemData.h"
#include "HelperFunctions.h"


std::vector<int> nodeVectorToIndexVector(const ProblemData& pd, const std::vector<MyGraph::Node>& solution);

// convert solution to vector with indices referring to destinations vector (containing the starting node as first element)
std::vector<int> nodeVectorToIndexVectorWithLimit(const ProblemData& pd, const std::vector<MyGraph::Node>& solution, std::size_t upperLimit);

std::vector<int> nodeVectorToIndexVectorWithStartNode(const ProblemData& pd, const std::vector<MyGraph::Node>& solution);

std::vector<int> nodeVectorToIndexVector(const ProblemData& pd, const std::vector<MyGraph::Node>& solution);

std::vector<MyGraph::Node> indexVectorToNodeVector(const ProblemData& pd, const int* indices,
                                                   const std::size_t numberOfElements,
                                                   const std::map<int, int>& indicesMap);

#endif
