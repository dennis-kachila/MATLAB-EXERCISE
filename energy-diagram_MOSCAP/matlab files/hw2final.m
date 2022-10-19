%{
         HW2
Developing simulation tool that would ask a user the type of the semiconductor (p or n type),
doping density, and surface potential in a silicon based MOS capacitor and then plot the energy diagram only in the semiconductor.
The code should ask the user to input the type of the semiconductor (p or n type), doping density N, and surface potential 𝜓s in a silicon based MOS capacitor and then plot the energy diagrams only in the semiconductor.
A code to solve Poisson’s equation for a MOS capacitor with a p-type semiconductor and n-type semiconductor.
The code should then plot the energy diagram only in the semiconductor.
Output: A figure (plot) with y-axis being Energy (eV), x-axis distance (μm) with x=0 being the surface of the semiconductor. 
The plot should have 4 curves for Ec, Ev, EF and Ei. 
%}

close all
clc

%{
The code should ask the user to input the type of the semiconductor (p or n type), doping density N, and surface potential 𝜓s in a silicon based MOS capacitor and then plot the energy diagrams only in the semiconductor.
%}

% Ask the user to input the type of the semiconductor (p or n type) if the user does not input p or n type, the default value is n type
type = input('Enter the type of the semiconductor (p or n type): ','s');
if type == 'p'
    type = 'p';
else
    type = 'n';
end

% Ask the user to input the doping density N, if no input is given, the default value is 1e16
N = input('Enter the doping density N: ');
if isempty(N)
    N = 1e16;
end

% Ask the user to input the surface potential 𝜓s, default value is 0.5
phis = input('Enter the surface potential 𝜓s: ');
if isempty(phis)
    phis = 0.5;
end


% set the constants
q = 1.6e-19; % charge of an electron (C)
k = 1.38e-23; % Boltzmann constant (J/K)
T = 300; % temperature (K)
ni = 1.5e10; % intrinsic carrier concentration (m^-3)
eps0 = 8.85e-12; % permittivity of free space (F/m)
epsSi = 11.7*eps0; % permittivity of silicon (F/m)
epsOx = 3.9*eps0; % permittivity of oxide (F/m)
tox = 10e-9; % thickness of oxide (m)
tsi = 200e-9; % thickness of silicon (m)
eSi = 11.7; % relative permittivity of silicon
eOx = 3.9; % relative permittivity of oxide
Eg = 1.12; % band gap of silicon (eV)
Eg0 = 1.17; % band gap of intrinsic silicon (eV)

% calculate the Fermi level
if type == 'p'
    EF = phis + (q*N*tsi)/(2*epsSi);
elseif type == 'n'
    EF = phis - (q*N*tsi)/(2*epsSi);
end

% calculate the intrinsic Fermi level using the Fermi-Dirac distribution
Ei = Eg0/2 + k*T/q*log(N/ni);


% calculate the conduction band using the Fermi-Dirac distribution
Ec = Eg/2 + k*T/q*log((N+ni*exp(-Eg/(2*k*T/q)))/(ni*exp(Eg/(2*k*T/q))));


% calculate the valence band using the Fermi-Dirac distribution
Ev = -Eg/2 + k*T/q*log((N+ni*exp(Eg/(2*k*T/q)))/(ni*exp(-Eg/(2*k*T/q))));


% solve Poisson’s equation for a MOS capacitor with a p-type semiconductor and n-type semiconductor
x = linspace(0,tsi,1000);
if type == 'p'
    y = (q*N*tsi)/(2*epsSi) - (q*N*x)/(epsSi);
elseif type == 'n'
    y = (q*N*tsi)/(2*epsSi) + (q*N*x)/(epsSi);
end

% plot the energy diagram only in the semiconductor
figure(1) % create a figure
plot(x,y,'k','LineWidth',2) %
hold on
plot([0 tsi],[Ec Ec],'r','LineWidth',2) % plot the conduction band
plot([0 tsi],[Ev Ev],'b','LineWidth',2) % plot the valence band
plot([0 tsi],[EF EF],'g','LineWidth',2) % plot the Fermi level
plot([0 tsi],[Ei Ei],'m','LineWidth',2) % plot the intrinsic Fermi level
hold off
xlabel('Distance (m)')
ylabel('Energy (eV)')
legend('Ei','Ec','Ev','EF')
title('Energy Diagram Only in the Semiconductor')
grid on



