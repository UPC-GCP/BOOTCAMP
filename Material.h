#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <vector>
#include <json/json.h>

struct MatPhys{
    double rho, lambda, cp, alpha;
};

class Material
{
private:

public:
    // Variables
    double T0{}, rho{}, lambda{}, cp{}, alpha{}, qV{};

    // Vectors
    std::vector<MatPhys> vMat{};

    // Constructor
    Material(double rho, double lambda, double cp, double source);
    // Material(Json::Value materials);
    
    // Functions
    void setInitialConditions(double initTemp);
    void setProperties(double rho, double lambda, double cp);
};

#endif