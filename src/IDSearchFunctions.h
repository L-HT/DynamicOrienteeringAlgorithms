#ifndef IDSEARCHFUNCTIONS_H
#define IDSEARCHFUNCTIONS_H

#include <Rcpp.h>

#include <lemon/list_graph.h>
#include <stdlib.h>
#include <lemon/lgf_writer.h>
#include <lemon/lgf_reader.h>
#include "ProblemData.h"

MyGraph::Node getNodeFromInternalID(const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const int& id);

#endif
