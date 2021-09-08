function [fB1956,Vp1,Vp2,Vs,rQp1,rQp2,rQs,GVp,GVs]=Biot_Gassmann_OneFluid(w,rhof,muf,Kf,rhos,Ks,phi,kappa,mu0,K0)
%% parameter
T=0.5*(1/phi+1);%Tortuosity
fB1956=log10(muf.*phi./kappa./rhof./T./2./pi);
%% 弹性模量中间变量
rhol=rhof.*T./phi+muf./(1i.*w.*kappa);%公式25
rho=(1-phi).*rhos+phi.*rhof;%岩石密度
rhoj=rho-rhof.^2./rhol;%计算速度用的密度 公式24
af=1-K0./Ks;%Biot系数 公式 4
M=Ks./(1-phi-K0./Ks+phi.*Ks./Kf);%耦合模量 公式 4
K=K0+M.*(af.^2);%饱和流体的岩石体积模量 公式 3
%% 计算原始Biot速度和耗散的计算 (Biot only)
vs=(mu0./rhoj).^(1/2);%横波复速度 
Vs=(real(1./vs)).^(-1);%横波实际速度 
%和纵波速度有关的量 (Biot only)
a1=(2.*af.*rhof-rho).*M-rhol.*(K+4./3.*mu0);% 公式 27
a0=(K0+4./3.*mu0).*M;
vp1=sqrt((-a1+(a1.^2-4.*rhoj.*rhol.*a0).^(1/2))./(2.*rhoj.*rhol));%复数速度，Vc上界限 
Vp1=(real(1./vp1)).^(-1);%快纵波
vp2=sqrt((-a1-(a1.^2-4.*rhoj.*rhol.*a0).^(1/2))./(2.*rhoj.*rhol));%复数速度，Vc下界限
Vp2=(real(1./vp2)).^(-1);%慢纵波
%求耗散 (Biot only)
rQp1=-imag(1./vp1).*2./real(1./vp1);
rQp2=-imag(1./vp2).*2./real(1./vp2);
rQs=-imag(1./vs).*2./real(1./vs);
%% 计算Gassmann速度
KGm=K0+(Ks./(1-phi-K0./Ks+phi.*Ks./Kf)).*((1-K0./Ks).^2);%Gassman饱和流体的岩石体积模量
GVp=((KGm+4/3*mu0)/rho)^0.5;
GVs=(mu0/rho)^0.5;
% GKu=rho.*(GVp.^2-4./3*GVs.^2);GNu=rho.*GVs.^2;%brine saturated rock
end