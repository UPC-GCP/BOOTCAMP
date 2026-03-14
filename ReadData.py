import os
import json
import decimal
import fractions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

###### Assignment Questions #####

def solverComparison(fileName:list, data:list):

    # Control
    sch, t, maxIter, res = [], [], [], []

    # File Loop
    for i, fName in enumerate(fileName):
        
        # Scheme
        iPos = fName.find('data_'); sch.append(fName[iPos+5:iPos+7])

        # Vectors
        t.append(data[i]["Time"]); maxIter.append(data[i]["MaxIter"]); res.append(data[i]["Residue"])

        # Clean
        t[i] = [float(x) for x in t[i][1:-1]]; maxIter[i] = [float(x) for x in maxIter[i][1:-1]]; res[i] = [float(x) for x in res[i][1:-1]]

    return (sch, t, maxIter, res)


def plotTempEvol(data:pd.DataFrame, x:float):
    
    # Plot temperature evolution at position "x"

    # Clean Vector
    iVec, xVec = [], []
    for i, xPos in enumerate(data.iloc[0]):
        if i == 0: continue
        iVec.append(i);
        try: 
            xVec.append(float(xPos))
        except:
            pass

    # Find Position
    xClean = [float(abs(fractions.Fraction(x) - fractions.Fraction(xVal))) for xVal in xVec[1:-1]]
    iFinal = np.argmin(xClean)
    
    # for xVal in xVec[1:-1]:
    #     xClean.append(float(abs(decimal.Decimal(x) - decimal.Decimal(xVal))))
    # xClean = [abs(float(decimal.Decimal(x)-decimal.Decimal(y))) for y in xVec[1:-1]]
    # iPos = int(np.where(xVec == x)[0])

    # Temperature Vector
    Temp = data["T" + str(iFinal)];
    Time = data["Time"]

    return [[float(x) for x in Time[1:-1]], Temp[1:-1]]


def plotTimeStep(data:pd.DataFrame, t:float):
    
    # Plot temperature profile at time "t" (s)

    # Clean Vector
    Time = data["Time"]; Time = [float(x) for x in Time[1:]]; Time = np.array(Time)
    
    # Find Position
    iPos = int(np.where(Time == t)[0])

    # Temperature Vector
    Temp = data.iloc[iPos+1][1:-3]
    xPos = data.iloc[0][1:-3]

    return (xPos, Temp)


def getTempVal(data:pd.DataFrame, x:float, t:float):
    
    # Find Time Step
    Time = data["Time"]; Time = [float(x) for x in Time[1:]]; Time = np.array(Time); iTime = int(np.where(Time == t)[0])
    
    # Find Position
    iVec, xVec = [], []
    for i, xPos in enumerate(data.iloc[0]):
        if i == 0: continue
        iVec.append(i); xVec.append(float(xPos))
    xVec = np.array(xVec); iPos = int(np.where(xVec == x)[0])

    # Value
    iVal = data.iloc[[iTime+1],[iPos+1]]
    
    print(f"Temperature @ x = {x:.3f} (m), t = {t:.2f} (s): ", iVal)


# Import Data
filename = '20260313075317_data_CG_crank-nicolson.csv'
data = pd.read_csv(filename)

def plotAnalytic(xVec:float, data:pd.DataFrame):

    Tl, Tr = 100, 20
    qV, V = 300000, 1
    lamb = 400

    C2 = Tl; C1 = 0.5 * qV * V / lamb + Tr - Tl
    
    yRet = [float(C2 + C1*x - 0.5 * qV * V * x * x / lamb) for x in xVec]

    return yRet


########## Plots ##########

# Plot A
xPos = 0.4
[xPlot, yPlot] = plotTempEvol(data, xPos)

fig = plt.figure(1); fig.set_figwidth(8); fig.set_figheight(4)
plt.plot(xPlot, yPlot, 'r')
plt.xscale('linear')
plt.grid(which='both', alpha=0.2); plt.title(f"Temperature evolution at position: {xPos:.3f} (m)" )
plt.xlabel('Time (s)'); plt.ylabel('Temperature (°C)')

# Plot B
tVal = 5000
[xPlot, yPlot] = plotTimeStep(data, tVal)
yAnal = plotAnalytic(xPlot, data)

fig = plt.figure(2); fig.set_figwidth(8); fig.set_figheight(4)
plt.plot(xPlot, yPlot, 'r', label='Model Solution')
plt.plot(xPlot, yAnal, 'b', label='Analytical Solution')
plt.grid(which='both', alpha=0.2); plt.title(f"Temperature profile at time: {tVal:.0f} (s)")
plt.xlabel('Position (m)'); plt.ylabel('Temperature (°C)')
plt.legend()


plt.show()
# Plot C
# getTempVal(data, 0.045, 20)



##### Mass Data Processing #####

def readMassResults():
    
    # TimeStep Study - What do I need to plot here?
    # dt vs Err (for each scheme, can separate it into arrays or can just get it into a single one a separate it afterwards)
    # N vs Err (for each scheme)

    # dt vs ComputationTime (for each scheme)
    # N vs ComputationTime (for each scheme)
    # CompTime not really needed

    # Mesh Study - What do I need to plot here?
    # How should I compare the meshes?

    # Control
    dir = os.getcwd(); strBase = "Case_";
    vName, vFile, vData = [], [], []

    # File List and .json Files
    for file in os.listdir(dir):
        if strBase in file and ".csv" in file:
            vName.append(str(file)); vFile.append(pd.read_csv(file, low_memory=False))
            with open(strBase + str(len(vFile)) + '.json') as srcData:
                vData.append(json.load(srcData)); srcData.close()

    # Control
    nt, nN, ii = 5, 4, 0
    xRet, yRet = [], []; xTemp, eTemp, iTemp, vTemp = [], [], [], []
    tSweep, NSweep = [], []; bN, bt = False, False

    # Processing Loop
    for i, name in enumerate(vName):
        
        # Control
        if i!=0 and i%9 == 0: ii+=1

        # Cases
        if i - ii*(nt + nN) < nt: # dtSweep
            
            # Check
            if bN:
                # Save Plots
                xRet.append(xTemp); yRet.append(eTemp)
                xRet.append(xTemp); yRet.append(iTemp)
                
                # Control
                bN = False
                xTemp, eTemp, iTemp = [], [], []
            
            # Populate Temp Arrays
            xTemp.append(vData[i]["timeStep"])

            vTemp = vFile[i]["Residue"][2:-1]
            vTemp = [float(x) for x in vTemp]
            eTemp.append(float(np.average(vTemp)))
            
            vTemp = vFile[i]["MaxIter"][2:-1]
            vTemp = [int(x) for x in vTemp]
            iTemp.append(float(np.average(vTemp)))

            # Control
            bt = True
            
        elif i - ii*(nt + nN) < nt+nN: # NSweep

            # Check
            if bt:
                # Save Plots
                xRet.append(xTemp); yRet.append(eTemp)
                xRet.append(xTemp); yRet.append(iTemp)

                # Control
                bt = False
                xTemp, eTemp, iTemp = [], [], []
            
            # Populate Temp Arrays
            xTemp.append(vData[i]["N"])

            vTemp = vFile[i]["Residue"][2:-1]
            vTemp = [float(x) for x in vTemp]
            eTemp.append(float(np.average(vTemp)))
            
            vTemp = vFile[i]["MaxIter"][2:-1]
            vTemp = [int(x) for x in vTemp]
            iTemp.append(float(np.average(vTemp)))
            
            # Control
            bN = True

        else:
            print(f"Didn't work: {i}")
        
    # Save to Array
    if bN:
        # Save Plots
        xRet.append(xTemp); yRet.append(eTemp)
        xRet.append(xTemp); yRet.append(iTemp)
        
        # Control
        bN = False
        xTemp, eTemp, iTemp = [], [], []
    
    print(len(xRet), len(yRet))

    # quit()

    return xRet, yRet

##### Mass Data #####
# [xPlot, yPlot] = readMassResults()
# 0, 1 = implicit (t) res, iter
# 2, 3 = implicit (N) res, iter
# 4, 5 = nicolson (t) res, iter
# 6, 7 = nicolson (N) res, iter
# 8, 9 = explicit (t) res, iter
# 10, 11 = explicit (N) res, iter
# Gets residue and iterations, need to include comparison with analytical solution as well

# fig = plt.figure(3); fig.set_figwidth(8); fig.set_figheight(4)
# plt.plot(xPlot[0], yPlot[0], label='implicit')
# plt.plot(xPlot[4], yPlot[4], label='crank-nicolson')
# plt.plot(xPlot[8], yPlot[8], label='explicit')
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# plt.xlabel(r'$\Delta t$'); plt.ylabel('Residue')

# fig = plt.figure(4); fig.set_figwidth(8); fig.set_figheight(4)
# plt.plot(xPlot[1], yPlot[1], label='implicit')
# plt.plot(xPlot[5], yPlot[5], label='crank-nicolson')
# plt.plot(xPlot[9], yPlot[9], label='explicit')
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# plt.xlabel(r'$\Delta t$'); plt.ylabel('# of Iterations')

# fig = plt.figure(5); fig.set_figwidth(8); fig.set_figheight(4)
# plt.plot(xPlot[2], yPlot[2], label='implicit')
# plt.plot(xPlot[6], yPlot[6], label='crank-nicolson')
# plt.plot(xPlot[10], yPlot[10], label='explicit')
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# plt.xlabel(r'N'); plt.ylabel('Residue')

# fig = plt.figure(6); fig.set_figwidth(8); fig.set_figheight(4)
# plt.plot(xPlot[3], yPlot[3], label='implicit')
# plt.plot(xPlot[7], yPlot[7], label='crank-nicolson')
# plt.plot(xPlot[11], yPlot[11], label='explicit')
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# plt.xlabel(r'N'); plt.ylabel('# of Iterations')


# # print(xPlot)
# # print(yPlot)

# plt.show()

# quit()


# ########## Data ##########

# # Import Data
# fileSolver, data = ['20260312043214_data_GS_crank-nicolson.csv', '20260312043225_data_CG_crank-nicolson.csv'], []; 

# for fName in fileSolver:
#     data.append(pd.read_csv(fName))


# ########## Plots ##########

# ##### Compare Linear Solvers #####
# [lPlot, xPlot, yIter, yRes] = solverComparison(fileSolver, data)
# fig = plt.figure(1); fig.set_figwidth(8); fig.set_figheight(4)
# for i, sch in enumerate(lPlot):
#     plt.plot(xPlot[i], yIter[i], label=sch)
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# # plt.title("Solver Comparison")
# plt.xlabel('Time (s)'); plt.ylabel('# Iterations')

# fig = plt.figure(2); fig.set_figwidth(8); fig.set_figheight(4)
# for i, sch in enumerate(lPlot):
#     plt.plot(xPlot[i], yRes[i], label=sch)
# plt.xscale('linear'); plt.yscale('linear')
# plt.grid(which='both', alpha=0.2); plt.legend()
# # plt.title("Solver Comparison")
# plt.xlabel('Time (s)'); plt.ylabel('Residue')








# End
# plt.show()