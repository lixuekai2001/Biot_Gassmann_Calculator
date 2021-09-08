clear all
clc
close all
%% input parameter
kappa=10^(-11);
rhos=2.65*10^3;
rhof=1*10^3;
phi=0.25;
K0=10*10^9; % dry frame bulk modulus
Ks=33*10^9;
Kf=2.2*10^9;
mu0=8.35*10^9;
mus=23*10^9;
muf=10^(-3); % shear viscosity
% xif=2.8*10^(-3); %bulk viscosity
% xif=0;
%% frequency
f=1:1:100;% 
ff(f)=log(10.^(-0.8+9./100.*f))./log(10);% 
w=2*pi*10.^(-0.8+9./100.*f);% 
%% calculation
% [v1,v2,v3,v4,q1,q2,q3,q4]=Muller_2011_function(w,rhof,muf,Kf,rhos,Ks,eta,kappa,mu0,K0,xif,mus);
[fB1956,Vp1,Vp2,Vs,rQp1,rQp2,rQs,GVp,GVs]=Biot_Gassmann_OneFluid(w,rhof,muf,Kf,rhos,Ks,eta,kappa,mu0,K0);
%% draw figures
figure(1)
plot(ff,Vp1,'r--','LineWidth',2)
hold on
plot_line(min(ff),GVp,max(ff),GVp,'b:',2)
set(gca,'linewidth',2)
xlabel('Frequency(Log(Hz))')
ylabel('Velocity (m/s)')
legend('Vp+£¨Biot£©','Vp+£¨Gassmann£©','location','best')
%
figure(2)
plot(ff,Vp2,'r--','LineWidth',2)
xlabel('Frequency(Log(Hz))')
ylabel('Velocity (m/s)')
legend('Vp-£¨Biot£©','location','best')
%
figure(3)
plot(ff,Vs,'r--','LineWidth',2)
hold on
plot_line(min(ff),GVs,max(ff),GVs,'b-',1)
xlabel('Frequency(Log(Hz))')
ylabel('Velocity (m/s)')
legend('Vs+£¨Biot£©','Vs+£¨Gassmann£©','location','best')
%
figure(5)
plot(ff,rQp1,'r--','LineWidth',2)
hold on
plot_line(fB1956,min(rQp1),fB1956,max(rQp1),'b:',1)
xlabel('Frequency(Log(Hz))')
ylabel('1/Q')
legend('1/Q (p+, Biot)','Biot frequency','location','best')
%
figure(6)
plot(ff,rQs,'r--','LineWidth',2)
hold on
plot_line(fB1956,min(rQs),fB1956,max(rQs),'b:',1)
xlabel('Frequency(Log(Hz))')
ylabel('1/Q')
legend('1/Q (s+, Biot)','Biot frequency','location','best')
%
figure(7)
plot(ff,rQp2,'r--','LineWidth',2)
hold on
plot_line(fB1956,min(rQp2),fB1956,max(rQp2),'b:',1)
xlabel('Frequency(Log(Hz))')
ylabel('1/Q')
legend('1/Q (p-, Biot)','Biot frequency','location','best')