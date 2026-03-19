// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"

Mesh::Mesh(int algo, double W, double H, double A, double xC, double kStr, double delta) {

    // Geometry
    this->W = W; this->H = H;

    // Mesh Parameters
    algorithm = algo;
    strength = A; centering = xC; kStrength = kStr; this->delta = delta;
    
}

bool Mesh::isFormula(std::string value){

    // Stringstream
    std::stringstream ss; ss << value;

    // Check
    float num = 0; ss >> num;

    // Return
    if (ss.good()){
        return true;
    } else if (num == 0 && value[0] != 0){
        return true;
    } else {
        return false;
    }

}

void Mesh::addBoundaryConditions(Json::Value boundaries){

    // Este módulo guarda los parámetros relevantes para cada BC para que después lo pueda usar el Discretizer.cpp y actualizar los boundary nodes.
    // Discretizer tiene que saber si puede solo usar el valor o si necesita actualizar el número. Puedo dejar el formato vector<double> por ahora y guardar un dato más ahí que me permita saber si hay que calcularlo.
    // OKOK ENTONCES ES ASÍ

    // 1. Resize boundaryExpr
    // 2. Read each BC and separate by type
        // Already done until here
    // 3. Mesh:
        // Current: boundaryConditions[i] = [type, position, value]
        // New: boundaryConditions[i] = [type, position, update, value]
        // IF isStr(boundaryConditions[i]["value"]) THEN update = 1 and save in boundaryExpr[i] ELSE update = 0
    // 4. Discretizer:
        // IF update = 1 THEN use boundaryConditions[i]["value"] = boundaryExpr[i](t)
        // 

    // Control
    boundaryExpr.resize(boundaries.size());

    // Boundary Conditions
    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++) {

        if (boundaries[i]["type"] == "Dirichlet") {

            if (isFormula(boundaries[i]["value"].asString())){

                // Save Data
                boundaryConditions.push_back({0, std::stod(boundaries[i]["position"].asString(), 0), 1, 0});
                boundaryExpr[i] = boundaries[i]["value"].asString();

            } else {

                // Save Data
                boundaryConditions.push_back({0, std::stod(boundaries[i]["position"].asString(), 0), 0, std::stod(boundaries[i]["value"].asString(), 0)});

            }
            

        } else if (boundaries[i]["type"] == "Neumann") {
            
            boundaryConditions.push_back({1, std::stod(boundaries[i]["position"].asString(), 0), std::stod(boundaries[i]["value"].asString(), 0), std::stod(boundaries[i]["side"].asString(), 0)});

        } else if (boundaries[i]["type"] == "Convection") {
            
            boundaryConditions.push_back({2, std::stod(boundaries[i]["position"].asString(), 0), std::stod(boundaries[i]["value"].asString(), 0), std::stod(boundaries[i]["side"].asString(), 0), std::stod(boundaries[i]["alpha"].asString(), 0)});

        } else {
            std::cerr << "Error: Invalid boundary condition type " << boundaries[i]["type"].asString() << std::endl;
        }

    }

}


void Mesh::calculateFaces(int cNode, int NSec, double x0, double x1) {

    // General
    double length = x1 - x0;

    // Face Positions
    if (algorithm == 0){
        // Face Positions 0: Bidirectional Non-uniform (A, xC)
        for (int i = cNode; i < cNode+NSec+1; i++) {
            xFaces[i] = x0 + (i-cNode) * length / NSec + strength * (centering - (i-cNode) * length / NSec) * (1 - (i-cNode)/NSec) * (i-cNode) / NSec;
        }
    } else if (algorithm == 1){
        // Face Positions 1: Unidirectional Non-uniform (Kappa)
        for (int i = cNode; i < cNode+NSec+1; i++) {
            xFaces[i] = x0 + pow(((i-cNode) * length / NSec), kStrength);
        }
    } else if(algorithm == 2){
        // Face Positions 2: Hyperbolic Tangent (Single Side)
        double A, B;
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta * ((static_cast<double>(i) - cNode) / NSec - 1)); B = tanh(delta);
            xFaces[i] = x0 + length * (1 + A / B);
        }
    } else if(algorithm == 3){
        // Face Positions 3: Hyperbolic Tangent (Double-Sided)
        double A, B;
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta*((static_cast<double>(i) - cNode)/NSec - 0.5));
            B = tanh(0.5 * delta);
            xFaces[i] = x0 + 0.5 * length * (1 + A/B);
        }
    }

}

void Mesh::generateMesh(Material& Mat, Json::Value sections) {
    
    // Control
    totNodes = 2;
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++) {
        totNodes += sections[i]["N"].asDouble();
    }
    xFaces.resize(totNodes-1); xNodes.resize(totNodes); xMat.resize(totNodes); qV.resize(totNodes);

    int cNode = 0;
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++) {

        // Material Vectors
        for (int j = 0; j < sections[i]["N"].asInt(); j++){
            xMat[1+j+cNode] = sections[i]["material"].asInt();
            qV[1+j+cNode] = sections[i]["qV"].asDouble();
        }

        // Face Nodes
        calculateFaces(cNode, sections[i]["N"].asInt(), sections[i]["x0"].asDouble(), sections[i]["x1"].asDouble());

        // Control
        cNode += sections[i]["N"].asInt();

    }

    // Boundaries
    xMat[0] = xMat[1]; xMat[xMat.size()-1] = xMat[xMat.size()-2];

    // CV Positions
    for (int i = 1; i < totNodes-1; i++) {
        xNodes[i] = 0.5 * (xFaces[i] + xFaces[i-1]);
    }
    xNodes.front() = xFaces.front(); xNodes.back() = xFaces.back();

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