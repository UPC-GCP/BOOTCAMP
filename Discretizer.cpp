// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"

Discretizer::Discretizer(std::string scheme, double endTime, double dt) {
    
    // Time Parameters
    this->scheme = scheme;
    this->endTime = endTime; this->dt = dt;

}

void Discretizer::setSchemeParameters(Material& Mat, Mesh& Msh){

    // Scheme Selection
    if (scheme == "explicit") {
        
        // Beta
        beta = 0;
        
        // Calculate Time-step
        std::vector<double> dtNew(Msh.totNodes, 0); double dtMin;
        for (int i = 0; i < Msh.N+2; i++){
            dtNew[i] = 0.5 * pow(Msh.deltaX[i], 2) / Mat.alpha;
        }

        // Update Time-step
        dtMin = *std::min_element(dtNew.begin()+1, dtNew.end()-1);
        if (dtMin < dt) {dt = dtMin;}
        
    } else if (scheme == "crank-nicolson") {

        // Beta
        beta = 0.5;

    } else if (scheme == "implicit") {

        // Beta
        beta = 1;

    } else {
        std::cerr << "Error: Invalid discretization scheme " << scheme << std::endl;
    }

}

void Discretizer::setBoundaryConditions(Material& Mat, Mesh& Msh){
    
    // Boundary Conditions
    int iPos;
    for (std::vector<double> bC : Msh.boundaryConditions) {

        // Control
        iPos = std::find(Msh.xNodes.begin(), Msh.xNodes.end(), bC[1]) - Msh.xNodes.begin();

        if (bC[0] == 0){
            
            // Dirichlet Coefficients
            Msh.TNodes[iPos] = bC[2];
            // Msh.ap[iPos] = 1; Msh.bp[iPos] = bC[2];

            // Control
            Msh.ignoreBC.push_back(iPos);

        } else if (bC[0] == 1) {
            
            // Neumann Coefficients
            if (!std::signbit(bC[3])) {
                Msh.matA[iPos].aw = - beta * Mat.lambda / Msh.dx[iPos-1];
                Msh.matA[iPos].ap = - Msh.matA[iPos].aw;
                Msh.bp[iPos] = bC[2] + (1 - beta) * (Mat.lambda * Msh.TNodes[iPos-1] / Msh.dx[iPos-1] - Mat.lambda * Msh.TNodes[iPos] / Msh.dx[iPos-1]);
            } else if (std::signbit(bC[3])){
                Msh.matA[iPos].ae = - beta * Mat.lambda / Msh.dx[iPos];
                Msh.matA[iPos].ap = - Msh.matA[iPos].ae;
                Msh.bp[iPos] = bC[2] + (1 - beta) * (Mat.lambda * Msh.TNodes[iPos+1] / Msh.dx[iPos] - Mat.lambda * Msh.TNodes[iPos] / Msh.dx[iPos]);
            } else {
                std::cerr << "Boundary side not specified correcly.\n";
            }

        } else if (bC[0] == 2) {

            // Convection Coefficients
            if (!std::signbit(bC[3])) {
                Msh.matA[iPos].aw = - beta * Mat.lambda / Msh.dx[iPos-1];
                Msh.matA[iPos].ap = - Msh.matA[iPos].aw + bC[4];
                Msh.bp[iPos] = bC[4] * bC[2] + (1 - beta) * (Mat.lambda * Msh.TNodes[iPos-1] / Msh.dx[iPos-1] - Mat.lambda * Msh.TNodes[iPos] / Msh.dx[iPos-1]);
            } else if (std::signbit(bC[3])){
                Msh.matA[iPos].ae = - beta * Mat.lambda / Msh.dx[iPos];
                Msh.matA[iPos].ap = - Msh.matA[iPos].ae + bC[4];
                Msh.bp[iPos] = bC[4] * bC[2] + (1 - beta) * (Mat.lambda * Msh.TNodes[iPos+1] / Msh.dx[iPos] - Mat.lambda * Msh.TNodes[iPos] / Msh.dx[iPos]);
            } else {
                std::cerr << "Boundary side not specified correcly.\n";
            }

        }

    }

}

// HASTA ACÁ ESTOY SEGURO DE QUE FUNCIONA BIEN, ESTA ES LA VERSIÓN CORRECTA DEL CÓDIGO PARA setBoundaryConditions

void Discretizer::setCoefficients(Material& Mat, Mesh& Msh){
    
    for (int i = 1; i < Msh.totNodes-1; i++) {

        // Coefficients A (Directo a matriz)
        Msh.matA[i].aw = -beta * Mat.lambda * Msh.Sw[i] / Msh.dx[i-1];
        Msh.matA[i].ae = -beta * Mat.lambda * Msh.Se[i] / Msh.dx[i];
        Msh.matA[i].ap = Mat.rho * Mat.cp * Msh.Vp[i] / dt - Msh.matA[i].aw - Msh.matA[i].ae;
        
        // Coefficients B
        Msh.bp[i] = Mat.qV * Msh.Vp[i] + Mat.rho * Mat.cp * Msh.Vp[i] * Msh.TNodes[i] / dt + (1 - beta) * Mat.lambda * (Msh.Sw[i]*Msh.TNodes[i-1]/Msh.dx[i-1] + Msh.Se[i]*Msh.TNodes[i+1]/Msh.dx[i] - (Msh.Sw[i]/Msh.dx[i-1] + Msh.Se[i]/Msh.dx[i])*Msh.TNodes[i]);

    }

}

void Discretizer::setRHS(Material& Mat, Mesh& Msh){

    for (int i = 1; i < Msh.totNodes-1; i++) {

        // Coefficients B
        Msh.bp[i] = Mat.qV * Msh.Vp[i] + Mat.rho * Mat.cp * Msh.Vp[i] * Msh.TNodes[i] / dt + (1 - beta) * Mat.lambda * (Msh.Sw[i]*Msh.TNodes[i-1]/Msh.dx[i-1] + Msh.Se[i]*Msh.TNodes[i+1]/Msh.dx[i] - (Msh.Sw[i]/Msh.dx[i-1] + Msh.Se[i]/Msh.dx[i])*Msh.TNodes[i]);
    
    }

}
