#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <json/json.h>

#include "Material.h"

struct Matrix{
    double ap=0, aw=0, ae=0;
};

// PENDIENTE
// Incluir un BC struct que me permita tener toda la información ordenada

class Mesh
{
private:

public:
    // Variables
    int totNodes{}, algorithm{};
    double W{}, H{}, strength{}, centering{}, kStrength{}, delta{};

    // Vectors
    std::vector<int> ignoreBC{}, xMat{};
    std::vector<std::vector<double>> boundaryConditions{};
    std::vector<std::string> boundaryExpr{};
    std::vector<double> xFaces{}, xNodes{}, TNodes{}, Sw{}, Se{}, dx{}, deltaX{}, Vp{}, bp{}, qV{};
    std::vector<Matrix> matA{};

    // Constructor
    Mesh(int algo, double W = 1, double H = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    bool isFormula(std::string value);
    void addBoundaryConditions(Json::Value boundaries);
    void calculateFaces(int cNode, int NSec, double x0, double x1);
    void generateMesh(Material& Mat, Json::Value sections);
};

#endif