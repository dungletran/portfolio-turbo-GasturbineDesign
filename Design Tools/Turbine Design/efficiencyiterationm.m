% run the designturbine.m function with initial isentropic efficiency, then
% the iteration will keep going until it doesn't change any more

clc
clearvars
eta_sai=0.88;
eta_sbi=0.88;
eta_sci=0.88;
eta_sdi=0.89;
eta_sei=0.7;
count=0;
while count<=1000
    designturbine(eta_sai,eta_sbi,eta_sci,eta_sdi,eta_sei);
    eta=xlsread('turbine data.xlsx',1,'B71:F71');
    if abs(eta(1)-eta_sai)>eps||abs(eta(2)-eta_sbi)>eps||abs(eta(3)-eta_sci)>eps||abs(eta(4)-eta_sdi)>eps||abs(eta(5)-eta_sei)>eps
        eta_sai=eta(1);
        eta_sbi=eta(2);
        eta_sci=eta(3);
        eta_sdi=eta(4);
        eta_sei=eta(5);
        count=count+1;
    else break
    end
end
    