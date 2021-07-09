// #ifndef BUDGETCHANGEHANDLING_H
// #define BUDGETCHANGEHANDLING_H
// 
// #include <Rcpp.h>
// 
// #include <iostream>
// #include <vector>
// #include <map>
// 
// #include "HelperFunctions.h"
// // #include "Solver.h"
// #include "ProblemData.h"
// #include "IDSearchFunctions.h"

/*
    * repair by memorization:
    *  -jedes mal bei evaluateSolution muss eine Tabelle gehalten werden -> neue Klasse?
    *      -TreeMap <Budget, solution>
    *      -lookup-Funktion: getBestSolution(Budget)
    *      -HashMap updaten wenn neue beste Lösung
    * 
    * no repair and restart:
    *      -als Solver neustarten? -> additionalLogData als neues Argument mitgeben?
    *      -nicht als Improvement-Heuristic, sondern als Solver?
    *          -Improvement-Heuristic würde viele Startlösungen haben 
    *              und sich ein bisschen mit "memorization" in die Quere kommen
    *      -
    *  repair by fixing:
    *      -nimm current/best und entferne Knoten: by maxConstraint-Violation, minTargetCriterion
    *    
    *  results muss man sich auch speichern, um Evaluationen zu vermeiden           
    */

// struct EvaluatedSolution {
//     std::vector<MyGraph::Node> solution_;
//     ResultData resultData_;
// };
// struct BudgetChangeHandler{
//     BudgetChangeHandler();
//     BudgetChangeHandler(int minBudget, int maxBudget, int tableSize);
//     int minBudget_;
//     int maxBudget_;
//     int tableSize_;
//     int typeOfHandling_;
//     
//     static const int HANDLING_NONE = 0;
//     static const int HANDLING_RESTART = 1;
//     static const int HANDLING_TABLE = 2;
//     static const int HANDLING_FIX_VIOLATION = 31;
//     static const int HANDLING_FIX_SPARINGLY = 32;
//     
//     std::vector<int> sampledInterval;
//     
//     std::map<int, EvaluatedSolution> bestSolutionsByBudget_;
//     
//     EvaluatedSolution getBestSolution(int budget);
//     
//     void initializeSampledTable(int minBudget, int maxBudget, int tableSize);
//     void updateTable(std::vector<MyGraph::Node> solution, ResultData resultData);
//     void fixInvalidSolutionViolation(ProblemData& problemData, 
//                                      std::vector<MyGraph::Node>& solution);
//     void fixInvalidSolutionSparingly(ProblemData& problemData, 
//                                      std::vector<MyGraph::Node>& solution);
//     ~BudgetChangeHandler();
// };


// #endif
                                            