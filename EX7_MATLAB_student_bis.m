close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);

path = "C:\Users\maeld\OneDrive\Bureau\EPFL\Physique_numerique\Exercice_7\"; %TODO insert path name
output= "test";
filename = path+output+"_f.out";
data_wave=load(filename);
filename = path+output+"_v.out";
velocity = load(filename);
filename = path+output+"_x.out";
data_x = load(filename);
time = data_wave(:,1);
wave = data_wave(:,2:end);

figure
pcolor(data_x,time,wave);shading interp;colorbar();xlabel("x [m]");ylabel("t [s]");

figure
contourf(data_x, time, wave);
xlabel('x [m]');
ylabel('t [s]');
