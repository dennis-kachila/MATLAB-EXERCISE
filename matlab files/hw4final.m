%{
    HW4
    Ask the user to enter the metal work function and calculate the surface potential 𝜓s.
    Ask user to enter total trapped charges in the oxide and calculate the Vfb.
    use the Vfb to develop code for plotting the energy diagram in the semiconductor and the oxide.

%}


close all
clc

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


% solve Poisson’s equation for a MOS capacitor with a p-type semiconductor and n-type semiconductor
x = linspace(0,tsi,1000);
if type == 'p'
    y = (q*N*tsi)/(2*epsSi) - (q*N*x)/(epsSi);
elseif type == 'n'
    y = (q*N*tsi)/(2*epsSi) + (q*N*x)/(epsSi);
end


%  user input to VG, thickness of the oxide d, and the oxide permittivity eOx else use default values
VG = input('Enter the gate voltage (V): ');
if isempty(VG)
    VG = 0.1;
end

d = input('Enter the thickness of the oxide d (m): ');
if isempty(d)
    d = 0.1e-6;
end

eOx = input('Enter the oxide permittivity (F/m): ');
if isempty(eOx)
    eOx = 3.9*8.85e-12;
end

% Ask the user to enter the metal work function and calculate the surface potential 𝜓s.
WF = input('Enter the metal work function (eV): ');
if isempty(WF)
    WF = 4.5;
end

% phis = WF - Eg/2;

% Ask user to enter total trapped charges in the oxide and calculate the Vfb.
Qox = input('Enter the total trapped charges in the oxide (C/m^2): ');
if isempty(Qox)
    Qox = 1e-10;
end
Vfb = (Qox*eOx)/(epsSi*d);
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

% use the Vfb to develop code for plotting the energy diagram in the semiconductor and the oxide.
figure(1) % create a figure
plot([0 tsi],[phis phis],'k','LineWidth',2) %
hold on
plot([0 tsi],[Ec Ec],'r','LineWidth',2) % plot the conduction band
plot([0 tsi],[Ev Ev],'b','LineWidth',2) % plot the valence band
plot([0 tsi],[EF EF],'g','LineWidth',2) % plot the Fermi level
plot([0 tsi],[Ei Ei],'m','LineWidth',2) % plot the instrinsic Fermi level
plot([0 tsi],[Vfb Vfb],'c','LineWidth',2) % plot the instrinsic Fermi level
hold off
xlabel('Distance (m)')
ylabel('Energy (eV)')
legend('Ei','Ec','Ev','EF')
title('Enegy Diagram using Vfb & Oxide Thickness')
