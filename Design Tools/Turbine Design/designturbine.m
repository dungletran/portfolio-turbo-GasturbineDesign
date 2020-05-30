% turbine design
% known : m, p_in, p_out, T_in, T_out, Dm_in, Dm_out, w angular speed.
% unkonwns: axial velocity, density, stage number, degree of reaction,
% stage load coefficient, flow coefficient. 
function []=designturbine(eta_sa,eta_sb,eta_sc,eta_sd,eta_se)

m=151.3; % kg/s combustion gas mass flow rate
p_in=873350; % Pa static pressure at first stage of turbine stator blade.
p_out=102200; % Pa static pressure at exit of turbine, to the atmosphere.
T_in=1222.7; % K temperature at first stage of turbine stator blade.
T_out=806.77; % K temperature at exit of turbine.
D_in=1.062; % m inlet mean diameter
D_out=1.08; % m exit mean diameter
w=469.35; % rad/s angular velocity
k=1.333; 
R=287;
cp=1148; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   preliminary estimation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_out_s=T_in*(p_out/p_in)^((k-1)/k); % isentropic exit temp
%eta=0.89; % total isentropic efficiency (T_in-T_out)/(T_in-T_out_s)
roh1a=p_in/R/T_in; % inlet density
roh_e=p_out/R/T_out; % exit density
% a b c d e corresponding stage 1, 2, 3, 4, 5
Da=D_in;

Db=(D_out-D_in)/4+D_in;

Dc=(D_out-D_in)/2+D_in;

Dd=(D_out-D_in)*3/4+D_in;

De=D_out;

% assume inlet velocity is axial, thus, alpha1 = 90. 
U3a=w*Da/2; % rotational speed for first stage rotor m/s assume U1=U2=U3 for the first stage of turbine
U3b=w*Db/2;
U3c=w*Dc/2;
U3d=w*Dd/2;
U3e=w*De/2; % rotational speed for last stage rotor m/s 
% last stage degree of reaction is 0
%r_exit=0;
%r=0.5; % for first,second and third fourth stage degree of reaction
%assume phi = 0.7 as the stage flow coefficient, all the other parameters
%can be solved.
% assume first four stage load coefficient lambdaa
lambda1=1.7;
r=0.5;
phi=0.7;
alpha1a=90; % first stage stator inlet angle
alpha2=acot((1-r+lambda1/2)/phi)*180/pi; % first four stages rotor inlet angle
alpha3=acot((1-r-lambda1/2)/phi)*180/pi+180; % first four stages rotor outlet angle
alpha1b=alpha3; % stage 2,3,4 stator inlet angle
beta2=180-alpha3; % first 4 stages rotor metal angle
beta3=180-alpha2;% first 4 stages rotor metal angle
Vaxa=U3a*phi; 
Vaxb=U3b*phi;
Vaxc=U3c*phi;
Vaxd=U3d*phi;
Vaxe=U3e*phi;
V3a=Vaxa/sin(alpha3*pi/180);
V3b=Vaxb/sin(alpha3*pi/180);
V3c=Vaxc/sin(alpha3*pi/180);
V3d=Vaxd/sin(alpha3*pi/180);

V1a=Vaxa;
V1b=V3a;
V1c=V3b;
V1d=V3c;
V1e=V3d;

W2a=V3a;
W2b=V3b;
W2c=V3c;
W2d=V3d;

U2a=U3a;
U2b=U3b;
U2c=U3c;
U2d=U3d;

V2a=Vaxa/sin(alpha2*pi/180);
V2b=Vaxb/sin(alpha2*pi/180);
V2c=Vaxc/sin(alpha2*pi/180);
V2d=Vaxd/sin(alpha2*pi/180);

W3a=V2a;
W3b=V2b;
W3c=V2c;
W3d=V2d;
% solve inlet area
A1a=m/roh1a/Vaxa;
h_b1a=A1a/pi/Da;
% velocity triangle solved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% as lambda=lm/U^2=cp*(T01-T03)/U^2
T1a=T_in;
T01a=T1a+V1a^2/2/cp;
H01a=cp*T01a;
T2a=T01a-V2a^2/2/cp;
T02a=T01a;
M2a=V2a/sqrt(k*R*T2a);

T03a=T01a-lambda1*U3a^2/cp;
T3a=T03a-V3a^2/2/cp;
H03a=cp*T03a;
H02a=H01a;
p1a=p_in;
p01a=p1a*(T01a/T1a)^(k/(k-1));
%T2as=T2a-0.05*V2a^2/2/cp;

eta_n=0.89; % define the nozzle isentropic efficiency
T2as=(1-(1-T2a/T01a)/eta_n)*T01a;
p2a=p01a/((T01a/T2as)^(k/(k-1)));
roh2a=p2a/R/T2a;
p02a=p2a+V2a^2*roh2a/2;
A2a=m/roh2a/Vaxa;
h_b2a=A2a/pi/Da;
% the isentropic efficiency at first stage is eta, which is also the whole thermal
% efficiency
T3as=T1a-(T1a-T3a)/eta_sa;
p3a=p1a*(T3as/T1a)^(k/(k-1));
roh3a=p3a/R/T3a;
p03a=p3a+V3a^2/2*roh3a;
A3a=m/roh3a/Vaxa;
h_b3a=A3a/pi/Da;
dT0a=T01a-T03a;
%%% first stage design finished  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for second stage the velocity triangle should be the same.the only
% difference is the absolute value since U change and the thermodynamic properties.


p1b=p3a;
roh1b=roh3a;
A1b=m/roh1b/Vaxb;
h_b1b=A1b/pi/Db;
T1b=T3a;
T01b=T1b+V1b^2/2/cp;
H01b=cp*T01b;
T2b=T01b-V2b^2/2/cp;
T02b=T01b;
H02b=H01b;
M2b=V2b/sqrt(k*R*T2b);
T2bs=(1-(1-T2b/T01b)/eta_n)*T01b;
T03b=T01b-lambda1*U3b^2/cp;
T3b=T03b-V3b^2/2/cp;
H03b=cp*T03b;
p01b=p1b*(T01b/T1b)^(k/(k-1));
p2b=p01b/((T01b/T2bs)^(k/(k-1)));
roh2b=p2b/R/T2b;
p02b=p2b+V2b^2*roh2b/2;
A2b=m/roh2b/Vaxb;
h_b2b=A2b/pi/Db;
% the isentropic efficiency at second stage is eta, which is also the whole thermal
% efficiency
T3bs=T1b-(T1b-T3b)/eta_sb;
p3b=p1b*(T3bs/T1b)^(k/(k-1));
roh3b=p3b/R/T3b;
p03b=p3b+V3b^2/2*roh3b;
A3b=m/roh3b/V3b;
h_b3b=A3b/pi/Db;
dT0b=T01b-T03b;

% for third stage
p1c=p3b;
roh1c=roh3b;
A1c=m/roh1c/Vaxc;
h_b1c=A1c/pi/Dc;
T1c=T3b;
T01c=T1c+V1c^2/2/cp;
H01c=cp*T01c;
T2c=T01c-V2c^2/2/cp;
T02c=T01c;
H02c=H01c;
M2c=V2c/sqrt(k*R*T2c);
T2cs=(1-(1-T2c/T01c)/eta_n)*T01c;
T03c=T01c-lambda1*U3c^2/cp;
T3c=T03c-V3c^2/2/cp;
H03c=cp*T03c;
p01c=p1c*(T01c/T1c)^(k/(k-1));
p2c=p01c/((T01c/T2cs)^(k/(k-1)));
roh2c=p2c/R/T2c;
p02c=p2c+V2c^2*roh2c/2;
A2c=m/roh2c/Vaxc;
h_b2c=A2c/pi/Dc;
% the isentropic efficiency at second stage is eta, which is also the whole thermal
% efficiency
T3cs=T1c-(T1c-T3c)/eta_sc;
p3c=p1c*(T3cs/T1c)^(k/(k-1));
roh3c=p3c/R/T3c;
p03c=p3c+V3c^2/2*roh3c;
A3c=m/roh3c/Vaxc;
h_b3c=A3c/pi/Dc;
dT0c=T01c-T03c;

%%%% for fourth stage %%%%%%%%%%
p1d=p3c;
roh1d=roh3c;
A1d=m/roh1d/Vaxd;
h_b1d=A1d/pi/Dd;
T1d=T3c;
T01d=T1d+V1d^2/2/cp;
H01d=cp*T01d;
T2d=T01d-V2d^2/2/cp;
T02d=T01d;
H02d=H01d;
M2d=V2d/sqrt(k*R*T2d);
T2ds=(1-(1-T2d/T01d)/eta_n)*T01d;
T03d=T01d-lambda1*U3d^2/cp;
T3d=T03d-V3d^2/2/cp;
H03d=cp*T03d;
p01d=p1d*(T01d/T1d)^(k/(k-1));
p2d=p01d/((T01d/T2ds)^(k/(k-1)));
roh2d=p2d/R/T2d;
p02d=p2d+V2d^2*roh2d/2;
A2d=m/roh2d/Vaxd;
h_b2d=A2d/pi/Dd;
% the isentropic efficiency at second stage is eta, which is also assumed as the whole thermal
% efficiency
T3ds=T1d-(T1d-T3d)/eta_sd;
p3d=p1d*(T3ds/T1d)^(k/(k-1));
roh3d=p3d/R/T3d;
p03d=p3d+V3d^2/2*roh3d;
A3d=m/roh3d/Vaxd;
h_b3d=A3d/pi/Dd;
dT0d=T01d-T03d;


%%%%%% last stage %%%%%%%%
re=0.6;
roh1e=roh3d;
% based on the above lamda_exit=1
lambdae=2*(1-re);
alpha1e=alpha3;
alpha2e=acot((1-re+lambdae/2)/phi)*180/pi;
beta2e=acot((lambdae/2-re)/phi)*180/pi+180;
alpha3e=90;
beta3e=acot(-(lambdae/2+re)/phi)*180/pi+180;
V3e=Vaxe/sin(alpha3e*pi/180);
% Vaxe=Vaxd;
% V3e=Vaxd;
W2e=Vaxe/sin(beta2e*pi/180);
U2e=U3e;
V2e=Vaxe/sin(alpha2e*pi/180);
W3e=sqrt(Vaxe^2+U3e^2);

% solve inlet area
A1e=m/roh1e/Vaxe;
h_b1e=A1e/pi/De;
% velocity triangle solved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% as lambda=lm/U^2=cp*(T01-T03)/U^2
T1e=T3d;
T01e=T1e+V1e^2/2/cp;
H01e=cp*T01e;
T2e=T01e-V2e^2/2/cp;
T02e=T01e;
H02e=H01e;
M2e=V2e/sqrt(k*R*T2e);
T2es=(1-(1-T2e/T01e)/0.8)*T01e;
T03e=T01e-lambdae*U3e^2/cp;
T3e=T03e-V3e^2/2/cp;
H03e=cp*T03e;
p1e=p3d;
p01e=p1e*(T01e/T1e)^(k/(k-1));
p2e=p01e/((T01e/T2es)^(k/(k-1)));
roh2e=p2e/R/T2e;
p02e=p2e+V2e^2*roh2e/2;
A2e=m/roh2e/Vaxe;
h_b2e=A2e/pi/De;
% the isentropic efficiency at first stage is eta, which is also the whole thermal
% efficiency
T3es=T1e-(T1e-T3e)/eta_se;
p3e=p1e*(T3es/T1e)^(k/(k-1));
roh3e=p3e/R/T3e;
p03e=p3e+V3e^2/2*roh3e;
A3e=m/roh3e/Vaxe;
h_b3e=A3e/pi/De;
dT0e=T01e-T03e;
% output result
p=[p1a p2a p3a p1b p2b p3b p1c p2c p3c p1d p2d p3d p1e p2e p3e];
T=[T1a T2a T3a T1b T2b T3b T1c T2c T3c T1d T2d T3d T1e T2e T3e];
roh=[roh1a roh2a roh3a roh1b roh2b roh3b roh1c roh2c roh3c roh1d roh2d roh3d roh1e roh2e roh3e];
p0=[p01a p02a p03a p01b p02b p03b p01c p02c p03c p01d p02d p03d p01e p02e p03e];
T0=[T01a T02a T03a T01b T02b T03b T01c T02c T03c T01d T02d T03d T01e T02e T03e];
H0=[H01a H02a H03a H01b H02b H03b H01c H02c H03c H01d H02d H03d H01e H02e H03e];
thermo=[roh;T;p;p0;T0];
D=[Da Db Dc Dd De];
V=[V1a V2a V3a V1b V2b V3b V1c V2c V3c V1d V2d V3d V1e V2e V3e];
W=[0 W2a W3a 0 W2b W3b 0 W2c W3c 0 W2d W3d 0 W2e W3e];
U=[0 U2a U3a 0 U2b U3b 0 U2c U3c 0 U2d U3d 0 U2e U3e];
alpha=[alpha1a alpha2 alpha3 alpha1b alpha2 alpha3 alpha1b alpha2 alpha3 alpha1b alpha2 alpha3 alpha1b alpha2e  alpha3e];
beta=[0 beta2 beta3 0 beta2 beta3 0 beta2 beta3 0 beta2 beta3 0 beta2e  beta3e];
v=[V; W; U;alpha;beta];
h_b=[h_b1a h_b2a h_b3a h_b1b h_b2b h_b3b h_b1c h_b2c h_b3c h_b1d h_b2d h_b3d h_b1e h_b2e h_b3e];
A=[A1a A2a A3a A1b A2b A3b A1c A2c A3c A1d A2d A3d A1e A2e A3e];
dT0=[dT0a dT0b dT0c dT0d dT0e];
for i=1:5
    Dt(3*(i-1)+1:i*3)=h_b(3*(i-1)+1:i*3)+D(i);
    Dh(3*(i-1)+1:i*3)=-h_b(3*(i-1)+1:i*3)+D(i);
end
geometry=[Dh; Dt; A; h_b];
designp=[r r r r re;lambda1 lambda1 lambda1 lambda1 lambdae];
% save it to excel
xlswrite('turbine data.xlsx',v,1,'B3:P7');
xlswrite('turbine data.xlsx',geometry,1,'B8:P11');
xlswrite('turbine data.xlsx',thermo,2,'B3:P7');
xlswrite('turbine data.xlsx',D,3,'B4:F4');
xlswrite('turbine data.xlsx',designp,3,'B2:F3');


