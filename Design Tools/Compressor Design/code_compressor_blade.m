%                     MEEN-646: Module 2                                     
%           Design of Subsonic Compressor Blade M2-C                                                        


% Objective:                                                                              
%    Develop a design software that enables you to generate subsonic compressor blades 
%    using NACA Camberline Equation                                                                                    
% Given Parameters                                                                        
%    - Base profile                
%    - Blade chord C
%    - Lift coefficient
%                                                                                         
% Instruction:                                                                            
%    Input: prin, K, H, t, s, CL, iZone selection

clc;
clear all; 


ln = @(x)(log(x));

%input:
iZone=1; % input 1, 2, 3 (1 for thinnest profile and 3 for thickest profile) 
prin=0; % Enter 1 if you want to print the output
% Enter your chord length data and lift coefficient data in K and H
K=[0.021470881].*1000; % chord length
H=-[-0.697243515]; % lift coefficient
% Enter value of t to choose ith value from matrix K and H
t=1;
ch = K(t); %chord length (created by Jita)
CL = H(t); %lift coefficient

%input parameters
C = 1; %blade chord (created by Tran)
s = 1; %spacing
%n_iter = 10000; %number of iteration

n_iter = 10000/10;

n_iter_b1 = 190/10;
n_iter_b2 = 6650/10;
n_iter_b3 = 9955/10; %for zone3 only
x_cam1=n_iter_b1/n_iter*C;
x_cam2=n_iter_b2/n_iter*C;
x_cam3=n_iter_b3/n_iter*C;%for zone3 only

%NACA Camberline Equation 
for i=1:1:n_iter
    x_cam(i) = i/n_iter*C;
    y_cam(i) = -C*CL/(4*pi)*((1-i/n_iter)*ln(1-i/n_iter)+i/n_iter*ln(i/n_iter));
    %camber line tangent angle
    v_cam(i)=atan(CL/(4*pi)*ln((1-i/n_iter)/(i/n_iter)));
    
    %for the 2nd camberline
    y_cam1(i) = y_cam(i)+s;
end

%blade thickness
if iZone == 1
    for i=1:1:n_iter_b1
        x(i)=i/n_iter;
        t(i)=C*(0.3419*x(i)^0.4929);%zone1
    end
    
    for i=(n_iter_b1+1):1:n_iter_b2
        x(i)=i/n_iter;
        t(i)=C*(-15.631*x(i)^6 + 38.563*x(i)^5 - 38.22*x(i)^4 + 19.934*x(i)^3 - 6.2802*x(i)^2 + 1.1333*x(i) + 0.0307);%zone1
    end
    
    for i=(n_iter_b2+1):1:n_iter
        x(i)=i/n_iter;
        t(i)=C*(75.656*x(i)^6 - 375.15*x(i)^5 + 774.1*x(i)^4 - 850.22*x(i)^3 + 524.07*x(i)^2 - 172.08*x(i) + 23.628);%zone1
    end
    
elseif iZone == 2
    for i=1:1:n_iter_b1
        x(i)=i/n_iter;
        t(i)=(C*(0.6128*x(i)^0.4937));%zone2
    end
    
    for i=(n_iter_b1+1):1:n_iter_b2
        x(i)=i/n_iter;
        t(i)=C*(-35.559*x(i)^6 + 83.97*x(i)^5 - 79.529*x(i)^4 + 39.519*x(i)^3 - 11.876*x(i)^2 + 2.0934*x(i) + 0.0531);%zone2
    end
    
    for i=(n_iter_b2+1):1:n_iter
        x(i)=i/n_iter;
        t(i)=C*(93.702*x(i)^6 - 455.68*x(i)^5 + 921.82*x(i)^4 - 991.96*x(i)^3 + 598.52*x(i)^2 - 192.34*x(i) + 25.931);%zone2 
    end
    
    if t(i)<0
       t(i)=0;
    end    
    
elseif iZone == 3
    for i=1:1:n_iter_b1
        x(i)=i/n_iter;
        t(i)=C*(0.8232*x(i)^0.4941);%zone3
    end
    
    for i=(n_iter_b1+1):1:n_iter_b2
        x(i)=i/n_iter;
        t(i)=C*(-56.476*x(i)^6 + 129.12*x(i)^5 - 118.24*x(i)^4 + 56.666*x(i)^3 - 16.456*x(i)^2 + 2.8703*x(i) + 0.0696);%zone3
    end
    
    for i=(n_iter_b2+1):1:n_iter_b3
        x(i)=i/n_iter;
        t(i)=C*(65.209*x(i)^6 - 309.36*x(i)^5 + 610.82*x(i)^4 - 641.17*x(i)^3 + 376.88*x(i)^2 - 118.1*x(i) + 15.713);%zone3
    end
    
    for i=(n_iter_b2+1):1:n_iter
        x(i)=i/n_iter;
        t(i)=C*(7.1679*x(i)^2 - 14.367*x(i) + 7.1992);%zone3
    end
    
    if t(i)<0
       t(i)=0;
    end 
    
end 
    

%suction side coordinate
for i=1:1:n_iter
    x_S(i) = x_cam(i) - (t(i)/2)*sin(v_cam(i));
    y_S(i) = y_cam(i) + (t(i)/2)*cos(v_cam(i));
    y_S1(i) = y_S(i) + s;
end

%pressure side coordinate
for i=1:1:n_iter
    x_P(i) = x_cam(i) + (t(i)/2)*sin(v_cam(i));
    y_P(i) = y_cam(i) - (t(i)/2)*cos(v_cam(i));
    y_P1(i) = y_P(i) + s;
end

% Scaling done to match the blade to actual chord length
x_cam=x_cam.*ch;y_cam=y_cam.*ch;
x_S=x_S.*ch; y_S=y_S.*ch;
x_P=x_P.*ch; y_P=y_P.*ch;
y_S1=y_S1.*ch;y_P1=y_P1.*ch;


% plot the blade
if prin==1
figure(1);
plot(x_cam,y_cam,'g')
hold on
plot(x_S,y_S,'r')
hold on
plot(x_P,y_P,'b')
hold on
plot(x_cam,y_cam1,'g')
hold on
plot(x_S,y_S1,'r')
hold on
plot(x_P,y_P1,'b')
hold on
xlim([-.1*ch 1.1*ch])
% ylim([-0.2 0.8])
ylim([-0.5*ch 0.5*ch])
end 

% Z matrix created to organize data for solidworks
Z(:,1)=x_S; Z(:,2)=y_S; Z(:,4)=x_P; Z(:,5)=y_P; Z(:,6)=0;

