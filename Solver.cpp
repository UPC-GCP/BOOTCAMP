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

Solver::Solver(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol){
    
    // Data
    maxIter = maxIterations; this->tolNum = tolNum; this->tolTemp = tolTime;

    // Timestamp
    time_t timeStamp = std::time(nullptr);
    struct tm datetime = *localtime(&timeStamp);
    
    // File Name
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    int iPos = fName.find(".json");
    fileName = oName + fName.substr(0, iPos) + "_" + fSol + "_" + scheme + ".csv";
    // PENDING: Move saved file to new location

}

void Solver::openFile(int N, std::vector<double> vNode, std::vector<double> vTemp, std::ofstream& file){

    // Control
    if (!file.is_open()){
        std::cerr << "Failed to open file. \n";
    }

    // File Headers
    file << "Time,";
    for (int i = 0; i < N; i++){
        file << "T" << i << ",";
    }
    file << "MaxIter,Residue,\n";

    file << "xNodes,";
    for (int i = 0; i < N; i++){
        file << vNode[i] << ",";
    }
    file << "MaxIter,Residue,\n";

    file << 0 << ",";
    for (int i = 0; i < N; i++){
        file << vTemp[i] << ",";
    }
    file << "0,0,\n";

}

void Solver::saveFile(int N, double t, std::vector<double> vTemp, std::ofstream& file){

    file << t << ",";
        for (int i = 0; i < N; i++) {
            file << vTemp[i] << ",";
        }
        file << lastIter << "," << lastRes << "\n";

}

void Solver::solve(std::vector<Matrix> matA, std::vector<double>& x, std::vector<double> matB, std::vector<int> ignoreBC){
    std::cout << "Virtual function \n";
}

void Solver::printNotes(std::ofstream& file, double tTot){

    // Info dump at the end of the file
    file << tTot << "\n";

}