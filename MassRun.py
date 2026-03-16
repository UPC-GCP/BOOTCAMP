import numpy as np
import json
import os
import subprocess

def createJsonFiles():

    # Count
    i = 0

    # Base Dictionary
    data = {
        "length": 1, "width": 1, "height": 1,
        "rho": 8960, "lambda": 400, "cp": 380,
        "T0": 30, "qV": 300000,
        "maxIterations": 2000, "tolNumeric": 1e-8,
        "endTime": 5000, "tolTemporal": 1e-3,
        "meshAlgorithm": 3, "strength": 0, "centering": 0.5, "kappa": 1, "delta": 0.001,
        "solver": "CG",

        "boundaries": [
            {"type": "Dirichlet", "position": "0", "value": "100"},
            {"type": "Dirichlet", "position": "1", "value": "20"},
            {"type": "Neumann", "position": "0", "value": "10000", "side": "-1"},
            {"type": "Neumann", "position": "1", "value": "-15000", "side": "1"},
            {"type": "Convection", "position": "0", "value": "80", "side": "-1", "alpha": "2000"},
            {"type": "Convection", "position": "1", "value": "10", "side": "1", "alpha": "4000"}
        ],

        "scheme": 'implicit', "N":50 , "timeStep": 0.5
    }

    # Control
    fBase = "Case_"

    # Ranges
    sSweep = ['implicit', 'crank-nicolson', 'explicit'] # "scheme"
    NSweep, N0 = [50, 100, 200, 400], 100 # "N"
    tSweep, t0 = [0.1, 0.5, 1, 2, 5], 0.5 # "timeStep"

    # Scheme Loop
    for scheme in sSweep:

        data.update({"scheme": scheme}); data.update({"N": N0}); data.update({"timeStep": t0})

        # Timestep Loop
        for dt in tSweep:

            # Update Value
            data.update({"timeStep": dt})

            # Print .json
            i = i+1; sName = fBase + str(i) + '.json'
            with open(sName, 'w') as file:
                json.dump(data, file, indent=4)
            print(f"Case {i}: {scheme} - N: {N0}, dt: {dt}")

        data.update({"timeStep": t0})

        # Nodes Loop
        for N in NSweep:

            # Update Value
            data.update({"N": N})

            # Print .json
            i = i+1; sName = fBase + str(i) + '.json'
            with open(sName, 'w') as file:
                json.dump(data, file, indent=4)
            print(f"Case {i}: {scheme} - N: {N}, dt: {t0}")

        data.update({"N": N0})

    print(f"{i} files created and saved to folder.")

def runSimulations():

    # Directory
    dName = os.path.dirname(os.path.realpath(__file__)) 

    # File Loop
    for i in range(1, 28):
        
        # Control
        print(f"Running Case {str(i)}:")
        sRun = 'MainSolver.exe Case_' + str(i) + '.json'
        cmd = os.path.join(dName, sRun)

        subprocess.run(cmd)

# Create .json Files
# createJsonFiles()

# Run Simulations
runSimulations()
    
