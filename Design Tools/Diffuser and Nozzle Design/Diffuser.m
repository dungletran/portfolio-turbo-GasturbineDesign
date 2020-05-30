clear all;

close all;



%Diffuser inlet ondition,Parameters from Turbine

p2 = 118300;                      %static pressure at inlet of diffuser

Cp = 1099;                        %heat capacity of air at T=804K

k = 1.354                         %heat capacity ratio

T2 = 804.87;

A2 = 1.6655;                      %annulus area at inlet of diffuser

rho2 = 0.51;

mfr = 151.3;                      %mass flow rate at the exit of turbine

R=287;                         

v2 = mfr/rho2/A2;

P20 = p2+0.5*rho2*v2*v2;

T20 = T2+v2^2/2/Cp;

c2 = sqrt(k*R*T2);

Ma2 = v2/c2;                      %Ma2 = 0.3124

p3 = 102200                       %diffuser exit static pressure

T3 = T20*(p3/P20)^((k-1)/k);



%calculation of nozzle inlet

D2_i = 0.958;                   

D3_i =0.958;                     

K = 2.5;                          % total pressure loss coefficient

P30 = P20-K*(P20-p2);             % exit total pressure

v3 = sqrt((P30-p3)*2/rho2);

rho3 = p3/R/T3;

c3 = sqrt(k*R*T3);

Ma3 = v3/c3;                      % exit Ma = 0.23

A3 = mfr/rho3/v3;

ARNd  = A3/A2;

D3_o = sqrt(4*A3/pi+D3_i^2);      %exit dia



