
close all
clc
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
Eg = 1.12; % band gap of silicon (eV)
Eg0 = 1.17; % band gap of intrinsic silicon (eV)

%{
The code ask the user to input the type of the semiconductor (p or n type), doping density N, and surface potential 𝜓s in a silicon based MOS capacitor and then plot the energy diagrams only in the semiconductor.
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

%{
            HW3
            Part 1
    limit the band bending to WDmax = 0.5V. If the band bending is greater than WDmax, the band bending is set to WDmax.
    Plot band bending for WDmax.
    Plot p(x ) and n(x ) as a function of x for positive x values only in the semiconductor.Use a semilog cordinate system for the plot.
    The code should take input to VG, thickness of the oxide d, and the oxide permittivity eOx then 
    calculate the surface potential 𝜓s and plot the energy diagram only in the semiconductor.
    
%}

%  limit the band bending to WDmax = 0.5V. If the band bending is greater than WDmax, the band bending is set to WDmax.
% set WDmax = 0.5V
WDmax = 0.5;
if phis > WDmax
    phis = WDmax;
end
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

% plot band bending for WDmax
figure(1) % create a figure
plot([0 tsi],[phis phis],'k','LineWidth',2) %
hold on
plot([0 tsi],[Ec Ec],'r','LineWidth',2) % plot the conduction band
plot([0 tsi],[Ev Ev],'b','LineWidth',2) % plot the valence band
plot([0 tsi],[EF EF],'g','LineWidth',2) % plot the Fermi level
plot([0 tsi],[Ei Ei],'m','LineWidth',2) % plot the instrinsic Fermi level
hold off
xlabel('Distance (m)')
ylabel('Energy (eV)')
legend('Ei','Ec','Ev','EF')
title('Energy Diagram Only in the Semiconductor using poison equation')
grid on

%{
    Plot p(x ) and n(x ) as a function of x for positive x values only in the semiconductor.
    Use a semilog cordinate system for the plot.
%}

%  Define the constants
q = 1.6e-19; %  electron charge
k = 1.38e-23; %  Boltzmann constant
T = 300; %  temperature in Kelvin
ni = 1.5e10; %  intrinsic carrier concentration in m^-3
Eg = 1.12; %  band gap in eV
Nc = 2.8e19; %  conduction band effective density of states in m^-3
Nv = 1.04e19; %  valence band effective density of states in m^-3
Ei = 0.026; %  intrinsic Fermi level in eV

%  Define the x values
x = 0:0.01:1;

%  Calculate the values of p and n
p = ni*exp((Eg/2 + Ei)*q/(k*T))*x;
n = ni*exp(-(Eg/2 - Ei)*q/(k*T))*x;

%  Plot the values of p and n
semilogy(x,p,x,n)
xlabel('x')
ylabel('p(x) and n(x)')
title('Plot of p(x) and n(x) as a function of x for positive x values only in the semiconductor')
legend('p(x)','n(x)')
grid on


%{
    Part 3
    The code should take user input to VG, thickness of the oxide d, and the oxide permittivity eOx then 
    calculate the surface potential 𝜓s and plot the energy diagram only in the semiconductor.   
%}

%  user input to VG, thickness of the oxide d, and the oxide permittivity eOx else use default values
VG = input('Enter the gate voltage (V): ');
if isempty(VG)
    VG = 0.1;
end
d = input('Enter the thickness of the oxide (m): ');
if isempty(d)
    d = 0.1e-6;
end
eOx = input('Enter the oxide permittivity (F/m): ');
if isempty(eOx)
    eOx = 3.9*8.85e-12;
end

% calculate the surface potential 𝜓s
phis = (q*VG)/(epsOx*d);

% calculate the intrinsic Fermi level using the surface potential 𝜓s
Ei = (q*phis)/(epsSi);

% calculate the conduction band
Ec = Ei + Eg/2;

% calculate the valence band
Ev = Ei - Eg/2;

% calculate the Fermi level
EF = (Ec + Ev)/2;

% plot the energy diagram only in the semiconductor
figure(3) % create a figure
plot([0 tsi],[phis phis],'k','LineWidth',2) %
hold on
plot([0 tsi],[Ec Ec],'r','LineWidth',2) % plot the conduction band
plot([0 tsi],[Ev Ev],'b','LineWidth',2) % plot the valence band
plot([0 tsi],[EF EF],'g','LineWidth',2) % plot the Fermi level
plot([0 tsi],[Ei Ei],'m','LineWidth',2) % plot the instrinsic Fermi level
hold off
xlabel('Distance (m)')
ylabel('Energy (eV)')
legend('Ei','Ec','Ev','EF')
title('Energy Diagram Only in the Semiconductor')
grid on


