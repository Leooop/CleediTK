clear
close all
% analytical solution for kinetic delay along a linear temperature path (T(x) = T0 + ga*x)
V = 3e-9; % slab velocity (m/s)
k0 = 1e-13; % reaction constant (s^-1)
Ea = 0; % activation energy J/mol
R = 8.314; % gas constant
T0 = 20; % temperature (K) at beginning of path
ga = 0.008; % temperature gradient (K/m) along the path

k0T = k0*exp(-Ea./(R*T0));

L = V/k0T; % characteristic length scale
EI = @(x) real(-expint(-x)); % exponential integral function

Nx = 1000;
x = linspace(0,10*L,Nx);
T = T0 + ga*x;
c = exp(-(k0/V)*(Ea/(ga*R))*(EI(-Ea./(R*T))  + (R*T/Ea).*exp(-Ea./(R*T))));
cT = exp(-(k0T/V)*x);

figure(1)
plot(x/L,c)
xlabel('distance normalized by (V/k0)')
ylabel('composition')
grid on 
hold on
plot(x/L,cT,'--')
