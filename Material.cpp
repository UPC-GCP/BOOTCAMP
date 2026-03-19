// Imports
#include <iostream>
#include <vector> 
#include <json/json.h>

// Self-Imports
#include "Material.h"

Material::Material(Json::Value materials){

    // if double - guardar directo en vMat[i].rho
    // if string - generar expresión y guardarla en varExpr, actualizar los valores antes de cada cálculo la función que actualice las propiedades también tiene que actualizar alpha aunque al final solo lo uso para el cálculo del dt en explicit. Tal vez no necesito tenerlo como array sino simplemente generarlo al momento de calcular el dt
    // dt se calcula usando propiedades termofísicas entonces si tengo propiedades variables y mi initial guess es muy lejano quiere decir que no voy a conseguir un dt funcional? Tipo todo va a faltar 
    
    // List
    vMat.resize(materials.size());

    // Materials
    for (Json::Value::ArrayIndex i = 0; i < materials.size(); i++){

        // Store Materials (Para Double)
        vMat[i].rho = materials[i]["rho"].asDouble();
        vMat[i].lambda = materials[i]["lambda"].asDouble();
        vMat[i].cp = materials[i]["cp"].asDouble();
        vMat[i].alpha = vMat[i].lambda / (vMat[i].rho * vMat[i].cp);

        // // Store Expression (Para String)
        // vMat[i].rhoExpr = materials[i]["rho"].asString();
        // vMat[i].lambdaExpr = materials[i]["lambda"].asString();
        // vMat[i].cpExpr = materials[i]["cp"].asString();

    }

}

void Material::setInitialConditions(double initTemp){

    // Initial Conditions
    T0 = initTemp;

}

