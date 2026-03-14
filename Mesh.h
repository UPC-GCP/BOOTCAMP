#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <json/json.h>

#include "Material.h"

struct Matrix{
    double ap=0, aw=0, ae=0;
};

class Mesh
{
private:

public:
    // Variables
    int N{}, totNodes{};
    double L{}, W{}, H{}, strength{}, centering{}, kStrength{}, delta{};

    // Vectors
    std::vector<int> ignoreBC{};
    std::vector<std::vector<double>> boundaryConditions{};
    std::vector<double> xFaces{}, xNodes{}, TNodes{}, Sw{}, Se{}, dx{}, deltaX{}, Vp{}, bp{};
    std::vector<Matrix> matA{};

    // Constructor
    Mesh(int N, double L, double W = 1, double H = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    void addBoundaryConditions(Json::Value boundaries);
    void generateMesh(Material& Mat, int algorithm);
};

#endif