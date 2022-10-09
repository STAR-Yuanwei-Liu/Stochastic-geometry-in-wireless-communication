%% SETTINGS %%
clear
clc
L = 0.75; % half length of RIS
fc = 1e9; % 1 GHz
c = 3*10^8; % speed of light
wavelength=3*10^8./fc;%wavelength
C_L=(wavelength/(4*pi))^2; %intercept of NLOS
BW=10*10^6; %bandwidth
Nf=10;%dB
Noise_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
Noise=10^((Noise_dbm-30)/10);
lambda=1/(300^2*pi); %density of PPP, the BSs
threshold=10^-2; %the threshold for SINR
threshold2=10^-2; %the threshold for SINR
r_C = 120; % connedted user distance to nearest user

Pt_dBm = -0:5:30; % transmit power dBm
P_t = 10.^((Pt_dBm-30)./10); % transmit power 
SNR = P_t./Noise; % transmit SNR 
SNR_dB = 10.*log10(SNR); % transmit SNR dB
num=1e5; % simulation times

R_BS = 4000; % the Range of BSs 
RL=25; % the radius of each cluster
Pout=zeros(num,1); % store space
m = 4; % gamma distribution
a_high = 0.6; % power allocation coefficient
a_low = 0.4; % power allocation coefficient
alpha = 4; % path loss coefficient 
alphat = 2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis 

rho_r = 0.4;
CP_t = zeros(1,length(SNR_dB));

for ii = 1 : length(SNR_dB)
    tic
    index = zeros(1,m);
    C_RIS_E = (L/4/pi)^2*(pi+sin(2*rho_r*pi)/(4*rho_r-12*rho_r^2+8*rho_r^3))/pi;

%     xishu = (alphat-2)*Noise/(pi*lambda*P_t(ii)*C_RIS_E);
%     fun = @(x,y) y./(x+xishu.*y.^(alphat).*x.^(alphat-1)).*exp(-pi*lambda.*x.^2);
%     threshold_E(ii) = 4/RL^2*(alphat-2)*integral2(fun,0,inf,0,RL);   
    threshold_E(ii) = 2*(alphat-2)*integral(@(x) exp(-pi*lambda.*x.^2)./x,0,inf);
%     threshold_E(ii) = (alphat-2)*(pi*lambda)^(-0.01/2)*gamma(0.01/2);

for n = 1:m
    miu_c = m*factorial(m).^(-1/m);
    miu1 = pi*lambda*(hypergeom([-2/alpha,m],1-2/alpha,-n*miu_c*threshold_E(ii)/m)-1);
    miu2 = (C_L*P_t(ii))^(-1)*n*miu_c*Noise*threshold_E(ii);
    miu4 =  n*miu_c*threshold*Noise/(a_high-a_low*threshold)/P_t(ii)/C_L;
    miu3 = pi*lambda*(hypergeom([-2/alpha,m],1-2/alpha,-n*miu_c*threshold/m/(a_high-a_low*threshold))-1);
%     index(n)=(-1)^(n+1)*nchoosek(m,n).*exp(-miu3.*r_C^2-miu4.*r_C^alpha);
    index(n)= -(-1)^(n+1)*nchoosek(m,n).*exp(-miu1.*r_C^2-miu2.*r_C^alpha)+(-1)^(n+1)*nchoosek(m,n).*exp(-miu3.*r_C^2-miu4.*r_C^alpha);

end
CP_c(ii) = sum(index);
toc
end

plot(SNR_dB,CP_c,'r-o');
hold on