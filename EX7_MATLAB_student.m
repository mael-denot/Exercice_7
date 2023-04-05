path = "/physiquenumerique/2022-2023/EX7/SOLUTION/"; %TODO insert path name 
filename = path+"test.out_f";
data_wave=load(filename);
filename = path+"test.out_v";
velocity = load(filename);
filename = path+"test.out_x";
data_x = load(filename);
time = data_wave(:,1);
wave = data_wave(:,2:end);

figure
pcolor(data_x,time,wave);shading interp;colorbar();xlabel("x [m]");ylabel("t [s]");
%% Egalement possible: contour ou contourf (tapez 'help contourf')

%% Modes propres, verification numerique 

A = 1; 
mode = 4;
f_analytic = A; % TODO: inserer la solution analytique du mode propre

figure;hold on;
plot(data_x,f_analytic)
plot(data_x,wave(1,1:end))
plot(data_x,wave(end,1:end))

%% Parameter scan: Excitation resonante
repertoire = "/physiquenumerique/2022-2023/EX7/SOLUTION/"; % TODO change as you need 
executable = './Exercice7.exe'; % TODO change as you need 
input = 'input'; % TODO change as you need 

omegamin=1; omegamax=10; nomega=5; % TODO: choose your own parameters
omega = linspace(omegamin, omegamax,nomega);
paramstr = 'omega'; 
param = omega;
output = cell(1, length(param));
Emax = zeros(1,length(param));
for i = 1:length(param)
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    cd /physiquenumerique/2022-2023/EX7/SOLUTION % TODO: choose your own path
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
    disp('Done.')
    filename = repertoire+output{i}+"_E";
    energy     = load(filename);
    Emax(i)    = max(max(energy));
end


figure
plot(omega,Emax) % TODO: xlabel, ylabel, etc...


%% compute the propagation velocity of the wave and its amplitude: TODO

%TODO: compare the WKB solutions with the numerical one

