import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation, matplotlib.ticker
import os, glob, argparse


def runSimulation(cppFile, params):
    '''
    Runs c++ simulation on input dictionary params and returns results as a matrix.
    Creates and then deletes two temporary files to pass parameters into and output data out of the c++ code.
    Results matrix contains output at each timestep as each row. Each variable is it's own column.

    Inputs:
        cppFile: String,    Name of c++ file
        params: dict,   Simulation Parameters

    Outputs:
        simulationData: np.array(nsteps+1, k),  matrix containing all output data across k variables
    '''

    # Set name of temporary output data file
    tempDataFileName = 'output_temp'
    params['output'] = tempDataFileName

    # Create temporary file of parameters to get passed into simulation. Deleted after use.
    tempParamsFileName = 'configuration_temp.in'
    with open(tempParamsFileName, 'w') as temp:
        for key in params.keys():
            temp.write("%s %s %s\n" % (key, "=", params[key]))

    # Run simulation from terminal, passing in parameters via temporary file
    os.system(' .\\' + cppFile + ".exe .\\" + tempParamsFileName)

    # Save output data in simData matrix
    simulationData_E   = np.loadtxt(tempDataFileName + '_E.out')
    simulationData_f   = np.loadtxt(tempDataFileName + '_f.out')
    simulationData_v = np.loadtxt(tempDataFileName + '_v.out')
    simulationData_x = np.loadtxt(tempDataFileName + '_x.out')
    time = simulationData_f[:,0]
    wave = simulationData_f[:,1:]

    # omega files (as of now, no omega files are outputted)

    # Delete temp files. Comment out if you want to keep them.
    os.remove(tempParamsFileName)
    os.remove(tempDataFileName + '_E.out')
    os.remove(tempDataFileName + '_f.out')
    # os.remove(tempDataFileName + '_v.out')
    # os.remove(tempDataFileName + '_x.out')

    return simulationData_E, simulationData_f, simulationData_v, simulationData_x, time, wave

programName = 'Exercice7_Denot_Wybaillie'

def plotSimulation(cppFile, params):
    E, f, v, x, time, psi = runSimulation(cppFile, params)

    # Ask for which variable to plot
    print("Which variable would you like to plot?")
    print("0: Exit")
    print("1: Energy")
    print("2: Wave")
    print("3: Velocity")
    print("4: All")
    choice = int(input("Enter your choice: "))
    
    # Plot the chosen variable
    if choice == 0:
        return
    elif choice == 1:
        plt.figure()
        plt.plot(E[:,0], E[:,1])
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.show()
    elif choice == 2:
        # do a contourf for f
        plt.figure()
        plt.contourf(x, time, psi)
        plt.xlabel("x")
        plt.ylabel("t")
        plt.show()
    elif choice == 3:
        plt.figure()
        plt.plot(x, v)
        plt.xlabel("x")
        plt.ylabel("v")
        plt.show()
    elif choice == 4:
        plt.figure()
        plt.plot(E[:,0], E[:,1])
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.show()
        plt.figure()
        plt.contourf(x, time, psi)
        plt.xlabel("x")
        plt.ylabel("t")
        plt.show()
        plt.figure()
        plt.plot(x, v)
        plt.xlabel("x")
        plt.ylabel("v")
        plt.show()
    else:
        print("Invalid choice. Try again.")
        plotSimulation(cppFile, params)

    # Ask if user wants to plot another variable
    print("Would you like to plot another variable?")
    print("0: No")
    print("1: Yes")
    choice = int(input("Enter your choice: "))
    if choice == 0:
        return
    elif choice == 1:
        plotSimulation(cppFile, params)
    else:
        print("Invalid choice. Try again.")
        plotSimulation(cppFile, params)

def contourplot(cppFile, params):
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")
    plt.show()

def compare_contourplot(cppFile, params):
    params['cb_gauche'] = 'harmonique'
    params['cb_droit'] = 'fixe'
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")

    params['cb_droit'] = 'libre'
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")

    plt.show()

#b) Limite de stabilite :
# Prendre un des cas de la partie precedente, avec un nombre d’intervalles nx donne. Verifier
# et illustrer que la solution devient instable des que |βCFL| > 1.

def stability(cppFile, params):
    params['cb_gauche'] = 'harmonique'
    params['cb_droit'] = 'fixe'
    params['CFL'] = 0.5
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")

    params['CFL'] = 1.0
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")

    params['CFL'] = 1.5
    E, f, v, x, time, psi = runSimulation(cppFile, params)
    plt.figure()
    plt.contourf(x, time, psi)
    plt.xlabel("x")
    plt.ylabel("t")

    plt.show()
    


parameters_basin = {
    'cb_gauche':'harmonique',
    'cb_droit':'fixe',
    'v_uniform' : 'true',
    'A':1.0,
    'omega':1.0,
    'tfin':4.0,

    'schema':'A',
    'Npoints':65,
    'minit':2,
    'initialization' : 'zero',
    'fmn':1.0, 

    'CFL':1.0,
    'output':'output.out',
    'n_stride':1,
    'ecrire_f':1,

    'hL' : 7.5e3,
    'hR' : 0.02e3,
    'h00': 2.0,
    'xa': 5e5,
    'xb': 10.5e5,
    'xL': 0.0,
    'xR': 12.0,
}


# plotting things
# parameters_basin['cb_droit'] = 'sortie'
# parameters_basin['cb_gauche'] = 'harmonique'
# parameters_basin['tfin'] = 2*12/1.0
# plotSimulation(programName, parameters_basin)

# 7.2a) réflexion aux bords
parameters_basin['omega'] = 7.5
parameters_basin['tfin'] = 1.5*12/3.5
# compare_contourplot(programName, parameters_basin)

parameters_basin['cb_gauche'] = 'harmonique'
parameters_basin['cb_droit'] = 'sortie'
# contourplot(programName, parameters_basin)

parameters_basin['cb_droit'] = 'fixe'
# contourplot(programName, parameters_basin)
# plotSimulation(programName, parameters_basin)

parameters_basin['cb_droit'] = 'fixe'
# contourplot(programName, parameters_basin)

# b) Limite de stabilite :

parameters_basin['cb_gauche'] = 'harmonique'
parameters_basin['cb_droit'] = 'fixe'
parameters_basin['CFL'] = 0.5
stability(programName, parameters_basin)