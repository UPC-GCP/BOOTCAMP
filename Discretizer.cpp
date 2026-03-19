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
#include "libFormulaParser.h"

Discretizer::Discretizer(std::string scheme, double endTime, double dt) {
    
    // Time Parameters
    this->scheme = scheme;
    this->endTime = endTime; this->dt = dt;

}

// double getExprValue(){

//     // This function will receive an expression for the corresponding equation and the control variable and must return the result 



// }

double Discretizer::calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX) {

    // Denominator
    double A = 0;
    for (int i = 0; i < lambda.size(); i++){
        A += (deltaX[i] / 2) / lambda[i];
    }
    
    return dPF / A;

}

void Discretizer::setSchemeParameters(Material& Mat, Mesh& Msh){

    // Scheme Selection
    if (scheme == "explicit") {
        
        // Beta
        beta = 0;
        
        // Calculate Time-step
        std::vector<double> dtNew(Msh.totNodes, 0); double dtMin;
        for (int i = 0; i < Msh.totNodes; i++){
            dtNew[i] = 0.5 * pow(Msh.deltaX[i], 2) / Mat.vMat[Msh.xMat[i]].alpha;
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

void Discretizer::setBoundaryConditions(Material& Mat, Mesh& Msh, double t){
    
    // Boundary Conditions
    int iPos;
    for (std::vector<double> bC : Msh.boundaryConditions) {

        // type, position, update, value

        // Control
        iPos = std::find(Msh.xNodes.begin(), Msh.xNodes.end(), bC[1]) - Msh.xNodes.begin(); double lamb;

        // Check Type
        if (bC[0] == 0){
            
            // Update Value
            if (bC[2] == 1){
                bC[3] = calcFormulaExpression(Msh.boundaryExpr[std::find(Msh.boundaryConditions.begin(), Msh.boundaryConditions.end(), bC) - Msh.boundaryConditions.begin()], t);
            }

            // Dirichlet Coefficients
            Msh.TNodes[iPos] = bC[3];
            // Msh.ap[iPos] = 1; Msh.bp[iPos] = bC[2];

            // Control
            Msh.ignoreBC.push_back(iPos);

        } else if (bC[0] == 1) {
            
            // Neumann Coefficients
            if (!std::signbit(bC[3])) {
                
                // Thermal Conductivity
                lamb = Mat.vMat[Msh.xMat[iPos]].lambda;

                // Coefficients
                Msh.matA[iPos].aw = - beta * lamb / Msh.dx[iPos-1];
                Msh.matA[iPos].ap = - Msh.matA[iPos].aw;
                Msh.bp[iPos] = bC[2] + (1 - beta) * (lamb * Msh.TNodes[iPos-1] / Msh.dx[iPos-1] - lamb * Msh.TNodes[iPos] / Msh.dx[iPos-1]);

            } else if (std::signbit(bC[3])){
                
                // Thermal Conductivity
                lamb = Mat.vMat[Msh.xMat[iPos]].lambda;

                // Coefficients
                Msh.matA[iPos].ae = - beta * lamb / Msh.dx[iPos];
                Msh.matA[iPos].ap = - Msh.matA[iPos].ae;
                Msh.bp[iPos] = bC[2] + (1 - beta) * (lamb * Msh.TNodes[iPos+1] / Msh.dx[iPos] - lamb * Msh.TNodes[iPos] / Msh.dx[iPos]);

            } else {
                std::cerr << "Boundary side not specified correcly.\n";
            }

        } else if (bC[0] == 2) {

            // Convection Coefficients
            if (!std::signbit(bC[3])) {

                // Thermal Conductivity
                lamb = Mat.vMat[Msh.xMat[iPos]].lambda;

                // Coefficients
                Msh.matA[iPos].aw = - beta * lamb / Msh.dx[iPos-1];
                Msh.matA[iPos].ap = - Msh.matA[iPos].aw + bC[4];
                Msh.bp[iPos] = bC[4] * bC[2] + (1 - beta) * (lamb * Msh.TNodes[iPos-1] / Msh.dx[iPos-1] - lamb * Msh.TNodes[iPos] / Msh.dx[iPos-1]);

            } else if (std::signbit(bC[3])){

                // Thermal Conductivity
                lamb = Mat.vMat[Msh.xMat[iPos]].lambda;

                // Coefficients
                Msh.matA[iPos].ae = - beta * lamb / Msh.dx[iPos];
                Msh.matA[iPos].ap = - Msh.matA[iPos].ae + bC[4];
                Msh.bp[iPos] = bC[4] * bC[2] + (1 - beta) * (lamb * Msh.TNodes[iPos+1] / Msh.dx[iPos] - lamb * Msh.TNodes[iPos] / Msh.dx[iPos]);

            } else {
                std::cerr << "Boundary side not specified correcly.\n";
            }

        }

    }

}

void Discretizer::setCoefficients(Material& Mat, Mesh& Msh){
    
    // Harmonic Mean
    double lambw, lambe;
    
    // Interior Nodes
    for (int i = 1; i < Msh.totNodes-1; i++) {

        // Harmonic Mean
        lambw = calcHarmonicMean(Msh.dx[i-1], {Mat.vMat[Msh.xMat[i-1]].lambda, Mat.vMat[Msh.xMat[i]].lambda}, {Msh.deltaX[i-1], Msh.deltaX[i]});
        lambe = calcHarmonicMean(Msh.dx[i], {Mat.vMat[Msh.xMat[i]].lambda, Mat.vMat[Msh.xMat[i+1]].lambda}, {Msh.deltaX[i], Msh.deltaX[i+1]});

        // Coefficients A (Directo a matriz)
        Msh.matA[i].aw = -beta * lambw * Msh.Sw[i] / Msh.dx[i-1];
        Msh.matA[i].ae = -beta * lambe * Msh.Se[i] / Msh.dx[i];
        Msh.matA[i].ap = Mat.vMat[Msh.xMat[i]].rho * Mat.vMat[Msh.xMat[i]].cp * Msh.Vp[i] / dt - Msh.matA[i].aw - Msh.matA[i].ae;
        
        // Coefficients B
        Msh.bp[i] = Msh.qV[i] * Msh.Vp[i] + Mat.vMat[Msh.xMat[i]].rho * Mat.vMat[Msh.xMat[i]].cp * Msh.Vp[i] * Msh.TNodes[i] / dt + (1 - beta) * (lambw*Msh.Sw[i]*Msh.TNodes[i-1]/Msh.dx[i-1] + lambe*Msh.Se[i]*Msh.TNodes[i+1]/Msh.dx[i] - (lambw*Msh.Sw[i]/Msh.dx[i-1] + lambe*Msh.Se[i]/Msh.dx[i])*Msh.TNodes[i]);

    }

}

void Discretizer::setRHS(Material& Mat, Mesh& Msh){

    // Harmonic Mean
    double lambw, lambe;

    // Interior Nodes
    for (int i = 1; i < Msh.totNodes-1; i++) {

        // Harmonic Mean
        lambw = calcHarmonicMean(Msh.dx[i-1], {Mat.vMat[Msh.xMat[i-1]].lambda, Mat.vMat[Msh.xMat[i]].lambda}, {Msh.deltaX[i-1], Msh.deltaX[i]});
        lambe = calcHarmonicMean(Msh.dx[i], {Mat.vMat[Msh.xMat[i]].lambda, Mat.vMat[Msh.xMat[i+1]].lambda}, {Msh.deltaX[i], Msh.deltaX[i+1]});

        // Coefficients B
        Msh.bp[i] = Msh.qV[i] * Msh.Vp[i] + Mat.vMat[Msh.xMat[i]].rho * Mat.vMat[Msh.xMat[i]].cp * Msh.Vp[i] * Msh.TNodes[i] / dt + (1 - beta) * (lambw*Msh.Sw[i]*Msh.TNodes[i-1]/Msh.dx[i-1] + lambe*Msh.Se[i]*Msh.TNodes[i+1]/Msh.dx[i] - (lambw*Msh.Sw[i]/Msh.dx[i-1] + lambe*Msh.Se[i]/Msh.dx[i])*Msh.TNodes[i]);
    
    }

}
