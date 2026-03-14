// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"

Mesh::Mesh(int N, double L, double W, double H, double A, double xC, double kStr, double delta) {

    // Geometry
    this->L = L; this->W = W; this->H = H;

    // Mesh Parameters
    this->N = N; strength = A; centering = xC; kStrength = kStr; this->delta = delta;
    
}

void Mesh::addBoundaryConditions(Json::Value boundaries){

    // Boundary Conditions
    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++) {

        if (boundaries[i]["type"] == "Dirichlet") {

            boundaryConditions.push_back({0, std::stod(boundaries[i]["position"].asString(), 0), std::stod(boundaries[i]["value"].asString(), 0)});

        } else if (boundaries[i]["type"] == "Neumann") {
            
            // PENDIENTE
            boundaryConditions.push_back({1, boundaries[i]["position"].asDouble(), boundaries[i]["value"].asDouble(), std::stod(boundaries[i]["side"].asString(), 0)});

        } else if (boundaries[i]["type"] == "Convection") {
            
            // PENDIENTE
            boundaryConditions.push_back({2, boundaries[i]["position"].asDouble(), boundaries[i]["value"].asDouble(), std::stod(boundaries[i]["side"].asString(), 0)});

        } else {
            std::cerr << "Error: Invalid boundary condition type " << boundaries[i]["type"].asString() << std::endl;
        }

    }

}

void Mesh::generateMesh(Material& Mat, int algo) {
    
    // Control
    totNodes = N + 2; 
    xFaces.resize(totNodes-1); xNodes.resize(totNodes);

    // Face Positions
    if (algo == 0){
        // Face Positions 0: Bidirectional Non-uniform (A, xC)
        for (int i = 0; i < totNodes-1; i++) {
            xFaces[i] = i * L / N + strength * (centering - i * L / N) * (1 - i/N) * i / N;
        }
    } else if (algo == 1){
        // Face Positions 1: Unidirectional Non-uniform (Kappa)
        for (int i = 0; i < totNodes-1; i++) {
            xFaces[i] = pow((i * L / N), kStrength);
        }
    } else if(algo == 2){
        // Face Positions 2: Hyperbolic Tangent (Single Side)
        double A, B;
        for (int i = 0; i < totNodes-1; i++){
            A = tanh(delta * (static_cast<double>(i) / N - 1)); B = tanh(delta);
            xFaces[i] = L * (1 + A / B);
        }
    } else if(algo == 3){
        // Face Positions 3: Hyperbolic Tangent (Double-Sided)
        double A, B;
        for (int i = 0; i < totNodes-1; i++){
            A = tanh(delta*(static_cast<double>(i)/N - 0.5));
            B = tanh(0.5 * delta);
            xFaces[i] = 0.5 * L * (1 + A/B);
        }
    }

    // CV Positions
    for (int i = 1; i < totNodes-1; i++) {
        xNodes[i] = 0.5 * (xFaces[i] + xFaces[i-1]);
    }
    xNodes.front() = 0; xNodes.back() = L;

    // Geometry
    Sw.resize(totNodes, W * H); Se.resize(totNodes, W * H); dx.resize(totNodes, 0); deltaX.resize(totNodes, 0); Vp.resize(totNodes, 0);

    // Delta
    for (int i = 0; i < totNodes - 1; i++) {
        dx[i] = xNodes[i+1] - xNodes[i];
        if (i > 0) {
            deltaX[i] = xFaces[i] - xFaces[i-1];
            Vp[i] = deltaX[i] * W * H;
        }
    }

    // Temperature
    TNodes.resize(totNodes, Mat.T0);

    // Coefficients
    matA.resize(totNodes); bp.resize(totNodes, 0);

}