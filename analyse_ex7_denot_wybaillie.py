import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation, matplotlib.ticker
import os, glob, argparse

path = "C:/Users/maeld/OneDrive/Bureau/EPFL/Physique_numerique/Exercice_7/" #TODO insert path name 
filename = path+"test_f.out";
data_wave=np.loadtxt(filename)
filename = path+"test_v.out";
velocity = np.loadtxt(filename)
filename = path+"test_x.out";
data_x = np.loadtxt(filename)
time = data_wave[:,0]
wave = data_wave[:,1:]

plt.figure()
plt.contourf(data_x,time,wave)
plt.colorbar()
plt.xlabel("x [m]")
plt.ylabel("t [s]")
plt.show()

# Modes propres, verification numerique 

A = 1; 
mode = 1; # TODO: inserer le mode number
f_analytic = A*np.zeros(data_x.size); # TODO: inserer la solution analytique du mode propre

plt.figure()
plt.plot(data_x,f_analytic)
plt.plot(data_x,wave[0,:])
plt.plot(data_x,wave[-1,:])
plt.show()

#%%

# Parameter scan: Excitation resonante
repertoire = "C:/Users/maeld/OneDrive/Bureau/EPFL/Physique_numerique/Exercice_7/"; # TODO change as you need 
executable = './Exercice7_Denot_Wybaillie.exe'; # TODO change as you need 
input = 'input'; # TODO change as you need 

omegamin=1; omegamax=10; nomega=5; # TODO: choose your own parameters
omega = np.linspace(omegamin, omegamax,nomega);
paramstr = 'omega'; 
param = omega;
#output = cell(1, length(param));
output = []
Emax = np.zeros([len(param)]);
for i in range(len(param)):
    output.append(paramstr+ '='+ str(param[i])+ '.out');
    os.system("cd {}".format(repertoire)) # TODO: choose your own path


    os.system("pwd")
    # cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    cmd = 'Exercice7.exe input omega='+str(param[i]) +' output='+paramstr+ '='+str(param[i])+'.out'
    os.system("./{}".format(cmd))

    filename = repertoire+output[i]+"_E";
    energy     = np.loadtxt(filename);
    Emax[i]    = np.max(np.max(energy));


plt.figure()
plt.plot(omega,Emax) # TODO: xlabel, ylabel, etc...

#%%
## compute the propagation velocity of the wave and its amplitude: TODO

#TODO: compare the WKB solutions with the numerical one




