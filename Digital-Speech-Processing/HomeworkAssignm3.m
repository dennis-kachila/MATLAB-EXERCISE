%{
    Question 5.12
     The transfer function of a digital formant resonator is of the form
     V_k(z) = (1 − 2*|z_k|*cos(θ_k) + |z_k|(^2))/(1 − 2*|zk|*cos(θ_k)*z^(-1) + |z_k|^2*z^(−2)) 
     where |zk| = e^(−σ_k*T) and θ_k = 2πf_k*T.

     Required:
     1. plot the location of the poles of the transfer function in the z-plane
     2. Plot the corresponding analogue poles in the s-plane
     3. To write the difference equation of the digital formant network with three multipliers.
    
%}

close all;
clc;
syms z;

% 1. plot the location of the poles of the transfer function in the z-plane and plot the corresponding analogue poles in the s-plane


% set the parameters
T = 1/8000;
f_k = 500;
sigma_k = 0.1;



% calculate the poles
z_k = exp(-sigma_k*T);
theta_k = 2*pi*f_k*T;

% set the z transfer function of the digital formant resonator
V_k = (1 - 2*abs(z_k).*cos(theta_k) + abs(z_k).^2)./(1 - 2*abs(z_k).*cos(theta_k).*z.^(-1) + abs(z_k).^2.*z.^(-2));


% denominator of the transfer function
denominator = 1 - 2*abs(z_k).*cos(theta_k).*z.^(-1) + abs(z_k).^2.*z.^(-2);
% extract the coefficients from the denominator
d = [1, -2*abs(z_k).*cos(theta_k), abs(z_k).^2];

% numerator of the transfer function
numerator = 1 - 2*abs(z_k).*cos(theta_k) + abs(z_k).^2;
% extract the coefficients from the numerator which is 0
n = [0, 0, 0];
% call the matlab built-in function to plot the poles in the z-plane
figure(1);
zplane(n, d);

% 2. plotting the corresponding analogue poles in the s-plane
figure(2);
pzmap(n, d);




%{
Question 5.17
The shape of the glottal pulse from the vocal cords can be approximated by the impulse
response of a second-order filter with system function:
G(z) = az^(−1)/(1 − a(z^(−1)))^(2) 

given  0 < a < 1.

(a) Plot the glottal pulse model impulse response, g[n], for a = 0.95 and for a = 0.8.
(b) Plot the corresponding log magnitude response, 20 log10 (|G(exp(jω)|, in dB versus ω (or
versus f) for the two values of a used in part (a) of this problem.
(c) The effect of lip radiation can be modeled by a single zero at z = 1. Repeat part (b)
with the inclusion of this extra zero.

%}

a = 0.95;
% The transfer function is given as follows
G_z = a*z^(-1)/(1 - a*z^(-1))^2;


% plotting the impulse response for a = 0.95
a = 0.95;
% calculate the impulse response

% denominator is expressed as
denominator = (1 - a*z^(-1))^2;
% expanding the denominator
d = [1 -19/10 361/400] ;
% numerator is expressed as
numerator = a*z^(-1) * z^2;

n = [19/20 0];

% call matlab function for discrete time impulse response
figure(3);
impz(n, d, 15);


% plotting the impulse response for
a = 0.8;
denominator = (1 - a*z^(-1))^2;
collect(denominator)
expand(denominator)
% expanding the denominator
d = [1 -8/5 16/25] ;
% numerator is expressed as
numerator = a*z^(-1);
n = [4/5 0];

% call matlab function for discrete time impulse response
figure(4);
impz(n, d, 15);




% % plotting the magnitude response for a = 0.95
%   defining the range of frequencies
w = 0:0.001:pi;
n = [19/20 0];
d = [1 -19/10 361/400] ;
%  calculating the magnitude response
% calling matlab freqz function
figure(5);
[H, w] = freqz(n, d, w);
% absolute value of H
mag = abs(H);
% plooting the magnitude response
plot(w/2, 20*log10(mag));
title('Magnitude Response for a = 0.95');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;

% % plotting the magnitude response for a = 0.8
%   defining the range of frequencies
w = 0:0.001:pi;
% expanding the denominator
d = [1 -8/5 16/25] ;
% numerator is expressed as
n = [4/5 0];
%  calculating the magnitude response
% calling matlab freqz function
figure(6);
[H, w] = freqz(n, d, w);
% absolute value of H
mag = abs(H);
% plooting the magnitude response
plot(w/2, 20*log10(mag));
title('Magnitude Response for a = 0.8');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;

w = 0:0.01:pi;

% adding an extra zero to the numerator of transfer function numerator
% since z = 1  for a = 0.95
new_n = [19/20 -1 0];
d = [1 -19/10 361/400] ;

%   calculating the magnitude response using the new numerator
%H = abs(19/20*exp(-j*w) + 361/400*exp(-2*j*w) - 19/20*exp(-j*w));
[H, w] = freqz(new_n, d, w);
%   plotting the magnitude response
figure(7);
plot(w, 20*log10(abs(H)));
title('Magnitude Response for a = 0.95 with extra zero');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;

% adding an extra zero to the numerator of transfer function numerator
% since z = 1  for a = 0.95
new_n = [19/20 -1 0];
d = [1 -19/10 361/400] ;

%   calculating the magnitude response using the new numerator
[ H, w ] = freqz(new_n, d, 1000);
%   plotting the magnitude response
figure(8);
plot(w, 20*log10(abs(H)));
title('Magnitude Response for a = 0.95 with extra zero');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;


% plotting the magnitude response for a = 0.8 with the extra zero
%   defining the range of frequencies
w = 0:0.01:pi;

% adding an extra zero to the numerator of transfer function numerator
% since z = 1  for a = 0.8
new_n = [4/5 -1 0];
d = [1 -8/5 16/25] ;

%   calculating the magnitude response using the new numerator
[ H, w ] = freqz(new_n, d, w);
%   plotting the magnitude response
figure(9);
plot(w, 20*log10(abs(H)));
title('Magnitude Response for a = 0.8 with extra zero');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;
