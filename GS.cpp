// Imports
#include <iostream>
#include <vector>
#include <string>
#include <json/json.h>
#include <cmath>
#include <math.h>
#include <numeric>
#include <ctime>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Solver.h"
#include "GS.h"
#include "libArithmetic.h"

void GS::solve(std::vector<Matrix> matA, std::vector<double>& x, std::vector<double> matB, std::vector<int> ignoreBC) {

    // Control
    double tempErr = 1, sum; size_t n = matB.size();
    std::vector<double> xOld, vRes(n), tempB(n);

    // Gauss-Seidel Loop
    for (int k = 0; k < maxIter; k++){

        // Nodes
        for (size_t i = 0; i < n; i++){
            
            // Control
            if (std::count(ignoreBC.begin(), ignoreBC.end(), i)){continue;}
            xOld = x;

            // Update Values
            if (i > 0 && i < n-1){
                x[i] = (-matA[i].aw * x[i-1] - matA[i].ae * x[i+1] + matB[i]) / matA[i].ap;
            } else if (i == 0){
                x[i] = (-matA[i].ae * x[i+1] + matB[i]) / matA[i].ap;
            } else if (i == n){
                x[i] = (-matA[i].aw * x[i-1] + matB[i]) / matA[i].ap;
            }

        }

        // Error
        vRes = operCombLinVec(operProdMatVec(matA, x), matB, 1, -1);
        tempErr = std::sqrt(operDotProd(vRes, vRes));
        if (tempErr < tolNum){lastIter = k; lastRes = tempErr; break;}

    }

}