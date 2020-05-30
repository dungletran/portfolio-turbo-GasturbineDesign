clear all; close all; clc;

%Nozzle inlet ondition,Parameters from compressor

mfr = 150;                            %mass flow rate of air

k = 1.4;

R = 287;

D1_o = 1.4506;

D1_i = 0.958;

D0_i =D1_i;

A1 = pi/4*(D1_o^2-D1_i^2);            %nozzle exit annulus area

p1 = 98610; 

T1 = 288.21;

rho = p1/R/T1;

v1 = mfr/rho/A1;

c1 = sqrt(k*R*T1);

Ma1 = v1/c1;
 
P1 = p1+0.5*rho*v1*v1;                %nozzle exit total pressure

%calculation of nozzle inlet

p0 = 101325;

Ploss = 0.05*P1;                      %5percent pressure loss

v0 = sqrt((P1-p0-Ploss)/0.5/rho);

A0 = mfr/rho/v0;                      %nozzle inlet area

T0 = p0/rho/R;

Ma0 = v0/c1;                          %nozzle inlet Ma = 0.197

ARN  = A0/A1;                         %area ratio

D0_o = sqrt(4*A0/pi+D0_i^2);          %nozzle inlet dia