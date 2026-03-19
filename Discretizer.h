#ifndef DISCRETIZER_H_
#define DISCRETIZER_H_

#include <string>

#include "Material.h"
#include "Mesh.h"

class Discretizer
{
private:

public:
    // Variables
    std::string scheme{};
    double beta{}, endTime{}, dt{};    

    // Constructor
    Discretizer(std::string scheme, double endTime, double dt);
    
    // Functions
    // double getExprValue(); // Receive double and return double. exprtk function called from Msh or Mat
    double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX);
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void setBoundaryConditions(Material& Mat, Mesh& Msh, double t = 0);
    void setCoefficients(Material& Mat, Mesh& Msh);
    void setRHS(Material& Mat, Mesh& Msh);
};

#endif