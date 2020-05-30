%                     MEEN-646: Module 2                                     
%           Design of Subsonic Turbine Blade M2-T                                                        

% Objective:                                                                              
%    Develop a design software that enables you to generate subsonic turbine blades 
%                                                                                        
% Given Parameters                                                                        
%    - Generate a family of profile (alpha1=90, alpha2=160), (alpha1=45, alpha2=160)              
%    - Blade chord C                                                                                     
%                                                                                         
% Instruction:                                                                            
%    Input: alpha1, alpha2 & iZone 

clc;
clearvars
iZone=1;
one=[110.69	139.46];
second=[156.63	154.3];
%input:
t=1;
alpha1 = one(t); %in degree
alpha2 = second(t); %in degree
%convert to rad
alpha1 = (alpha1/180)*pi;
alpha2 = (alpha2/180)*pi;

C = 40; %chord
s = 30; %spacing
n_iter = 1000; %number of iteration for camberline

n_iter_b1 = 19;
n_iter_b2 = 665;
n_iter_b3 = 995; %for zone3 only
x_cam1=n_iter_b1/n_iter*C;
x_cam2=n_iter_b2/n_iter*C;
x_cam3=n_iter_b3/n_iter*C;%for zone3 only

%cascade stagger angle
gamma = atan(sin(alpha2)/(-1/2*sin(alpha1-alpha2)/sin(alpha1)+cos(alpha2))); 
C_ax = C*sin(gamma);


%define camber line equation:
x_p0 = 0; y_p0 = 0;
x_p2 = C; y_p2 = 0;

%determine P1 coordinates by consider triangle P0P1P3
%a_1 = 1/3*C/sin(alpha1)*sin(gamma); %P1P0 length
%b_1 = 1/3*C/sin(alpha1)*sin(pi-alpha1-gamma); %P1P3 length
%c_1 = 1/3*C; %P0P3 length
%p = (a_1+b_1+c_1)/2;
%area = sqrt(p*(p-a_1)*(p-b_1)*(p-c_1));

%y_p1 = 2*area/(1/3*C);
%x_p1 = y_p1/tan(pi-alpha1-gamma);

%determine P1 coordinates, formula given in the book (equation 10.40)
phi1=pi/2-alpha1+gamma;
phi2=pi/2+alpha2-gamma;
y_p1 = C*(cot(phi1)/(1+cot(phi1)/cot(phi2)));
x_p1 = C*(1/(1+cot(phi1)/cot(phi2)));



%Bezier Curve
for i=1:1:n_iter
    zeta(i) = i/n_iter;
    x_cam(i) = (1-zeta(i))^2*x_p0 + 2*(1-zeta(i))*zeta(i)*x_p1 + zeta(i)^2*x_p2;
    y_cam(i) = (1-zeta(i))^2*y_p0 + 2*(1-zeta(i))*zeta(i)*y_p1 + zeta(i)^2*y_p2;
    
    %camber line tangent angle
    v_cam(i)=atan((-2*(1-zeta(i))*y_p0 + 2*(1-2*zeta(i))*y_p1 + 2*zeta(i)*y_p2)/(-2*(1-zeta(i))*x_p0 + 2*(1-2*zeta(i))*x_p1 + 2*zeta(i)*x_p2));
    
    %for the 2nd camberline
    y_cam1(i) = y_cam(i)+s;
end

%blade thickness
if iZone == 1
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.3419*x(i)^0.4929);%zone1
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-15.631*x(i)^6 + 38.563*x(i)^5 - 38.22*x(i)^4 + 19.934*x(i)^3 - 6.2802*x(i)^2 + 1.1333*x(i) + 0.0307);%zone1
        elseif (x_cam(i)>x_cam2)
            t(i)=C*(75.656*x(i)^6 - 375.15*x(i)^5 + 774.1*x(i)^4 - 850.22*x(i)^3 + 524.07*x(i)^2 - 172.08*x(i) + 23.628);%zone1
        end 
        %t(i)=t(i)/2;
    end
elseif iZone == 2
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.6128*x(i)^0.4937);%zone2
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-35.559*x(i)^6 + 83.97*x(i)^5 - 79.529*x(i)^4 + 39.519*x(i)^3 - 11.876*x(i)^2 + 2.0934*x(i) + 0.0531);%zone2
        elseif (x_cam(i)>x_cam2)
            t(i)=C*(93.702*x(i)^6 - 455.68*x(i)^5 + 921.82*x(i)^4 - 991.96*x(i)^3 + 598.52*x(i)^2 - 192.34*x(i) + 25.931);%zone2 
        end
        
        if t(i)<0
           t(i)=0;
        end
    end
elseif iZone == 3
    for i=1:1:n_iter
        x(i)=x_cam(i)/C;
        if (x_cam(i)<x_cam1)
            t(i)=C*(0.8232*x(i)^0.4941);%zone3
        elseif (x_cam(i)>x_cam1) && (x_cam(i)<x_cam2)
            t(i)=C*(-56.476*x(i)^6 + 129.12*x(i)^5 - 118.24*x(i)^4 + 56.666*x(i)^3 - 16.456*x(i)^2 + 2.8703*x(i) + 0.0696);%zone3
        elseif (x_cam(i)>x_cam2) && (x_cam(i)<x_cam3)
            t(i)=C*(65.209*x(i)^6 - 309.36*x(i)^5 + 610.82*x(i)^4 - 641.17*x(i)^3 + 376.88*x(i)^2 - 118.1*x(i) + 15.713);%zone3
        elseif (x_cam(i)>x_cam3)
            t(i)=C*(7.1679*x(i)^2 - 14.367*x(i) + 7.1992);%zone3
        end
        
        if t(i)<0
           t(i)=0;
        end
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
    %x_P_test(i) = x_cam(i) - (t(i)/2)*sin(v_cam(i));
    y_P(i) = y_cam(i) - (t(i)/2)*cos(v_cam(i));
    y_P1(i) = y_P(i) + s;
end


%plot the blade
figure(1);
plot(x_cam,y_cam,'g')
hold on
plot(x_S,y_S,'r')
hold on
plot(x_P,y_P,'b')
hold on
%plot(x_P_test,y_P,'y')
%hold on
%for i=1:1000:(n_iter)
%    th = 0:pi/50:2*pi;
%    xunit = t(i)/2 * cos(th) + x_cam(i);
%    yunit = t(i)/2 * sin(th) + y_cam(i);
%    h = plot(xunit, yunit,'r');
%    hold on
%end
xf=[x_S';x_P'];
yf=[y_S';y_P'];
zf=zeros(2000,1);
Z=[xf yf zf];
% plot(x_cam,y_cam1,'g')
% hold on
% plot(x_S,y_S1,'r')
% hold on
% plot(x_P,y_P1,'b')
% hold on
% 
% x_S=x_S'; x_P=x_P';
% y_P1=y_P1'; y_S1=y_S1';
% 
% %plot(x_P_test,y_P1,'y')
% %hold on
% %for i=1000:1000:(n_iter-1000)
% %    th = 0:pi/50:2*pi;
% %    xunit = t(i)/2 * cos(th) + x_cam(i);
% %    yunit = t(i)/2 * sin(th) + y_cam1(i);
% %    h = plot(xunit, yunit,'r');
% %    hold on
% %end
% 
% xlim([-4 44])
% %ylim([-0.2 1])
% ylim([-8 12])
% 
% gamma*180/pi
