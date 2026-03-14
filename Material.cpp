// Imports
#include <iostream>
#include <vector> 
#include <json/json.h>

// Self-Imports
#include "Material.h"

Material::Material(double density, double conductivity, double specificHeat, double source) {

    // Thermophysical Properties
    rho = density; lambda = conductivity; cp = specificHeat; qV = source;
    alpha = lambda / (rho * cp);

}

// Material::Material(Json::Value materials){

//     // List
//     vMat.resize(materials.size());

//     // Materials
//     for (Json::Value::ArrayIndex i = 0; i < materials.size(); i++){

//         // Store Materials
//         vMat[i].rho = materials[i]["rho"].asDouble();
//         vMat[i].lambda = materials[i]["lambda"].asDouble();
//         vMat[i].cp = materials[i]["cp"].asDouble();
//         vMat[i].alpha = vMat[i].lambda / (vMat[i].rho * vMat[i].cp);

//     }

// }

// Material::Material(Json::Value materials){
//     Make it store values in array so you can change materials for different parts.
// }

void Material::setInitialConditions(double initTemp){

    // Initial Conditions
    T0 = initTemp;

}

void Material::setProperties(double density, double conductivity, double specificHeat){
    
    // Thermophysical Properties
    rho = density; lambda = conductivity; cp = specificHeat;
    alpha = lambda / (rho * cp);

}