#ifndef IDSEARCHFUNCTIONS_H
#define IDSEARCHFUNCTIONS_H

#include <Rcpp.h>

#include <list_graph.h>
#include <stdlib.h>
#include <lgf_writer.h>
#include <lgf_reader.h>
#include "ProblemData.h"

MyGraph::Node getNodeFromInternalID(const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const int& id);

#endif
