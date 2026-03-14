#ifndef SOLVER_H_
#define SOLVER_H_

#include <string>
#include <vector>
#include <fstream>

#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"

class Solver
{
private:
    
public:
    // Variables
    std::string fileName{};
    double maxIter{}, tolNum{}, tolTemp{}, lastIter{}, lastRes{};

    // Constructor
    Solver(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol);

    // File Functions
    void openFile(int N, std::vector<double> vNode, std::vector<double> vTemp, std::ofstream& file);
    void saveFile(int N, double t, std::vector<double> vTemp, std::ofstream& file);

    // Inheritance Functions
    virtual void solve(std::vector<Matrix> matA, std::vector<double>& x, std::vector<double> matB, std::vector<int> ignoreBC = {}) = 0;
    void printNotes(std::ofstream& file, double tTot);
};

#endif