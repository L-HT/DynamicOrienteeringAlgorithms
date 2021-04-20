#include <vector>
#include <algorithm>
#include <random>
#include <lemon/list_graph.h>

#include "Solver.h"
#include "HelperFunctions.h"
#include "CalculateLength.h"
#include "ProblemData.h"

using namespace lemon;

std::vector<MyGraph::Node> apply2OptIteration(Solver& solver, const std::vector<MyGraph::Node>& solution, Rcpp::NumericMatrix* distanceMatrix){
    std::vector<MyGraph::Node> bestSolution(solution.begin(), solution.end());
    double shortestLength = solver.evaluateSolution(solver.problemData_, solution, solver.targetCriterion_, false, false, distanceMatrix).length_;

    // waitForInput("start 2opt", DEBUG_ENABLED);
    // printNodeIdsOfVector(solver.problemData_, solution);
    for (std::size_t lowerLimit = 0; lowerLimit < solution.size() - 1; lowerLimit++){
        for (std::size_t upperLimit = lowerLimit+1; upperLimit < solution.size(); upperLimit++){
            std::vector<MyGraph::Node> tempSolution(solution.size(), lemon::INVALID);
            for (std::size_t i = 0; i < tempSolution.size(); i++){
                if (i >= lowerLimit && i <= upperLimit){
                    tempSolution[i] = solution[upperLimit-(i - lowerLimit)];
                } else {
                    tempSolution[i] = solution[i];
                }
            }
            // printNodeIdsOfVector(solver.problemData_, tempSolution);
            ResultData res = solver.evaluateSolution(solver.problemData_, tempSolution, solver.targetCriterion_, false, false, distanceMatrix);
            if (res.length_ < shortestLength){
                shortestLength = res.length_;
                bestSolution.assign(tempSolution.begin(), tempSolution.end());
            }
        }
    }
    return bestSolution;
}
