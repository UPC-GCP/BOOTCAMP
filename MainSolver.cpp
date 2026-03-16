// Imports
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <json/json.h>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <thread>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Solver.h"
#include "GS.h"
#include "CG.h"

Json::Value getParsedData(std::string fileName){
    
    // Open File
    std::ifstream file(fileName, std::ifstream::binary);

    // Filter
    if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << fileName << std::endl;
    return 1;
    }

    // Parsing
    Json::Value data;
    Json::CharReaderBuilder readerBuilder;
    std::string errs;
    Json::parseFromStream(readerBuilder, file, &data, &errs);

    // Close File
    file.close();

    return data;

}

int main(int argc, char* argv[]){

    // Time
    auto t1 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Initializing model ... \n";
    
    ////////// Config File //////////
    std::cout << "Reading data ... \n";

    // Read Config
    Json::Value data = getParsedData(argv[1]);
    std::cout << "Data parsed successfully. \n";



    ////////// Model Implementation //////////

    ///// Material /////
    std::cout << "Initializing material ... \n";

    // Create Material
    Material Mat(data["rho"].asDouble(), data["lambda"].asDouble(), data["cp"].asDouble(), data["qV"].asDouble());
    std::cout << "Thermophysical properties set. \n";

    // Initial Conditions
    Mat.setInitialConditions(data["T0"].asDouble());
    std::cout << "Initial conditions set. \n";


    ///// Mesh /////
    std::cout << "Initializing mesh ... \n";

    // Create Mesh
    Mesh Msh(data["N"].asInt(), data["length"].asDouble(), data["width"].asDouble(), data["height"].asDouble(), data["strength"].asDouble(), data["centering"].asDouble(), data["kappa"].asDouble(), data["delta"].asDouble());
    std::cout << "Geometry set. \n";

    // Add Boundary Conditions
    Msh.addBoundaryConditions(data["boundaries"]);
    std::cout << Msh.boundaryConditions.size() << " boundary conditions added. \n";

    std::cout << "Generating mesh ... \n";

    // Generate Mesh
    Msh.generateMesh(Mat, data["meshAlgorithm"].asInt());
    std::cout << "Mesh created with " << Msh.totNodes << " nodes. \n";


    ///// Discretizer /////
    std::cout << "Initializing discretizer ... \n";

    // Create Discretizer
    Discretizer Dsc(data["scheme"].asString(), data["endTime"].asDouble(), data["timeStep"].asDouble());

    // Discretizer Calculations
    Dsc.setSchemeParameters(Mat, Msh); std::cout << "Temporal parameters set. \n";
    Dsc.setBoundaryConditions(Mat, Msh); std::cout << "Boundary conditions set. \n";
    Dsc.setCoefficients(Mat, Msh); std::cout << "Discretized coefficiets set. \n";

    
    ///// Solver /////
    std::cout << "Initializing solver ... \n";

    // Create Solver
    Solver* Sol = nullptr;
    if (data["solver"] == "CG"){
        Sol = new CG(Dsc.scheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
    } else if (data["solver"] == "GS"){
        Sol = new GS(Dsc.scheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
    } else {
        std::cerr << "Error: Invalid linear solver selected " << data["solver"].asString() << std::endl;
    }

    // Open File 
    std::ofstream file(Sol->fileName);
    Sol->openFile(Msh.totNodes, Msh.xNodes, Msh.TNodes, file);

    std::cout << "File created. \n";



    ////////// Temporal Loop //////////
    std::cout << "Processing ... \n";

    for (double t = Dsc.dt; t <= Dsc.endTime; t += Dsc.dt){

        // Update Coefficients
        Dsc.setBoundaryConditions(Mat, Msh);
        Dsc.setRHS(Mat, Msh);

        // Solver
        Sol->solve(Msh.matA, Msh.TNodes, Msh.bp, Msh.ignoreBC);

        // Write Content
        Sol->saveFile(Msh.totNodes, t, Msh.TNodes, file);

        std::cout << "\r" << double(100 * t / Dsc.endTime) << " %";

    }
    std::cout << "\n";

    // Time
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> msDoub = t2 - t1;

    std::cout << "Time elapsed: " << msDoub.count()/1000 << "\n";

    // Notes
    Sol->printNotes(file, msDoub.count()/1000);

    // Close File
    file.close();

    std::cout << "File saved: " << Sol->fileName << "\n";

}

