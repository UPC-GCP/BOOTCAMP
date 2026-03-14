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
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void setBoundaryConditions(Material& Mat, Mesh& Msh);
    void setCoefficients(Material& Mat, Mesh& Msh);
    void setRHS(Material& Mat, Mesh& Msh);
};

#endif