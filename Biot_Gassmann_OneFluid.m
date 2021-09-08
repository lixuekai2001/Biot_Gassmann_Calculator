function [fB1956,Vp1,Vp2,Vs,rQp1,rQp2,rQs,GVp,GVs]=Biot_Gassmann_OneFluid(w,rhof,muf,Kf,rhos,Ks,phi,kappa,mu0,K0)
%% parameter
T=0.5*(1/phi+1);%Tortuosity
fB1956=log10(muf.*phi./kappa./rhof./T./2./pi);
%% ����ģ���м����
rhol=rhof.*T./phi+muf./(1i.*w.*kappa);%��ʽ25
rho=(1-phi).*rhos+phi.*rhof;%��ʯ�ܶ�
rhoj=rho-rhof.^2./rhol;%�����ٶ��õ��ܶ� ��ʽ24
af=1-K0./Ks;%Biotϵ�� ��ʽ 4
M=Ks./(1-phi-K0./Ks+phi.*Ks./Kf);%���ģ�� ��ʽ 4
K=K0+M.*(af.^2);%�����������ʯ���ģ�� ��ʽ 3
%% ����ԭʼBiot�ٶȺͺ�ɢ�ļ��� (Biot only)
vs=(mu0./rhoj).^(1/2);%�Შ���ٶ� 
Vs=(real(1./vs)).^(-1);%�Შʵ���ٶ� 
%���ݲ��ٶ��йص��� (Biot only)
a1=(2.*af.*rhof-rho).*M-rhol.*(K+4./3.*mu0);% ��ʽ 27
a0=(K0+4./3.*mu0).*M;
vp1=sqrt((-a1+(a1.^2-4.*rhoj.*rhol.*a0).^(1/2))./(2.*rhoj.*rhol));%�����ٶȣ�Vc�Ͻ��� 
Vp1=(real(1./vp1)).^(-1);%���ݲ�
vp2=sqrt((-a1-(a1.^2-4.*rhoj.*rhol.*a0).^(1/2))./(2.*rhoj.*rhol));%�����ٶȣ�Vc�½���
Vp2=(real(1./vp2)).^(-1);%���ݲ�
%���ɢ (Biot only)
rQp1=-imag(1./vp1).*2./real(1./vp1);
rQp2=-imag(1./vp2).*2./real(1./vp2);
rQs=-imag(1./vs).*2./real(1./vs);
%% ����Gassmann�ٶ�
KGm=K0+(Ks./(1-phi-K0./Ks+phi.*Ks./Kf)).*((1-K0./Ks).^2);%Gassman�����������ʯ���ģ��
GVp=((KGm+4/3*mu0)/rho)^0.5;
GVs=(mu0/rho)^0.5;
% GKu=rho.*(GVp.^2-4./3*GVs.^2);GNu=rho.*GVs.^2;%brine saturated rock
end