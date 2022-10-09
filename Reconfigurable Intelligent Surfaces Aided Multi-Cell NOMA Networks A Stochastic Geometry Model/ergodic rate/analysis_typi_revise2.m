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

Pt_dBm = -0:30:30; % transmit power dBm
P_t = 10.^((Pt_dBm-30)./10); % transmit power 
SNR = P_t./Noise; % transmit SNR 
SNR_dB = 10.*log10(SNR); % transmit SNR dB
num=1e5; % simulation times

R_BS = 4000; % the Range of BSs 
RL=25; % the radius of each cluster
Pout=zeros(num,1); % store space
m = 1; % gamma distribution
a_high = 0.6; % power allocation coefficient
a_low = 0.4; % power allocation coefficient
alpha = 4; % path loss coefficient 
alphat = 2.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis 

rho_r = 0.4;
CP_t = zeros(1,length(SNR_dB));
for ii = 1 : length(SNR_dB)
    threshold_E = (alpha-2)*P_t(ii)*C_L*r_C^(-alpha)/ ( 2*pi*lambda*P_t(ii)*C_L*r_C^(2-alpha) + (alpha-2)*Noise );
    index = zeros(1,m);
    C_RIS_E = (L/4/pi)^2*(pi+sin(2*rho_r*pi)/(4*rho_r-12*rho_r^2+8*rho_r^3))/pi;
for n = 1:m
    tic
    miu_t = m*factorial(m).^(-1/m);
    A = [threshold2/(a_high-threshold2*a_low); threshold/a_low; threshold_E];
    A2 = [threshold2/(a_high-threshold2*a_low); threshold/a_low];
    Gamma = max(A2);
    beta0 = @(y) n*miu_t*Gamma*Noise/P_t(ii)/C_RIS_E.*y.^alphat;
    beta2 = pi*lambda*hypergeom([-2/alpha,m],1-2/alpha,-n*miu_t*Gamma/m);
    fun =@(x,y) 2*pi*lambda*(-1)^(n+1)*nchoosek(m,n).*x.*exp(-beta0(y).*x.^alphat).*exp(-beta2.*x.^2)...
        .*2.*y./RL^2;
    index1(n) = integral3(@(x,y,z) fun(x,y)./(1+z),0,1000,0,RL,0,a_low*Gamma);
    
    beta01 = @(y,z) n*miu_t.*z./a_low*Noise/P_t(ii)/C_RIS_E.*y.^alphat;
    beta21 = @(z) pi*lambda*hypergeom([-2/alpha,m],1-2/alpha,-n*miu_t.*z./a_low/m);
    fun =@(x,y,z) 2*pi*lambda*(-1)^(n+1)*nchoosek(m,n).*x.*exp(-beta01(y,z).*x.^alphat).*exp(-beta21(z).*x.^2)...
        .*2.*y./RL^2;
    index2(n) = integral3(@(x,y,z) fun(x,y,z)./(1+z),0,1000,0,RL,a_low*Gamma,1000);
    index(n) = index1(n)+index2(n) ;
    toc
end
CP_t(ii) =1/log(2)* sum(index);

end

plot(SNR_dB,CP_t,'r-o');
hold on
% CP_t = [0.0753230649775675,0.306893797603351,1.03325213715032,2.70839716622219];