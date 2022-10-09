%% SETTINGS %%
clear
clc
L = 2:5:22; % half length of RIS
fc = 1e9; % 1 GHz
c = 3*10^8; % speed of light
wavelength=3*10^8./fc;%wavelength
C_L=(wavelength/(4*pi))^2; %intercept of NLOS
BW=10*10^6; %bandwidth
Nf=10;%dB
Noise_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
Noise=10^((Noise_dbm-30)/10);
lambda=1/(600^2*pi); %density of PPP, the BSs
threshold=10^-2; %the threshold for SINR
threshold2=10^-2; %the threshold for SINR
r_C = 120; % connedted user distance to nearest user

Pt_dBm = 3; % transmit power dBm
P_t = 10.^((Pt_dBm-30)./10); % transmit power 
SNR = P_t./Noise; % transmit SNR 
SNR_dB = 10.*log10(SNR); % transmit SNR dB
num=1e5; % simulation times

R_BS = 4000; % the Range of BSs 
RL= 25; % the radius of each cluster
Pout=zeros(num,1); % store space
m = 4; % gamma distribution
a_high = 0.6; % power allocation coefficient
a_low = 0.4; % power allocation coefficient
alpha = 4; % path loss coefficient 
alphat = 2.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis 

rho_r = 0.4;
CP_t = zeros(1,length(L));

for ii = 1 : length(L)
    
    threshold_E = (alpha-2)*P_t*C_L*r_C^(-alpha)/ ( 2*pi*lambda*P_t*C_L*r_C^(2-alpha) + (alpha-2)*Noise );
    index_1 = zeros(1,m);
    index_2 = zeros(1,m);
    index_sum = zeros(1,m);
    C_RIS_E = (L(ii)/4/pi)^2*(pi+sin(2*rho_r*pi)/(4*rho_r-12*rho_r^2+8*rho_r^3))/pi;
    
for n = 1:m
    tic;
    miu_t = m*factorial(m).^(-1/m);
    xishu_1 = @(x,y,z) n*miu_t.*z.*Noise.*y.^alphat.*x.^alphat./(P_t*C_RIS_E*a_low);
    A = [threshold2/(a_high-threshold2*a_low); threshold/a_low; threshold_E];
    A2 = [threshold2/(a_high-threshold2*a_low); threshold_E];
    Gamma = max(A2);
    fun_index =@(y,z) quadgk(@(t) t-t.*(1+n*miu_t/a_low/m.*y.^alphat.*z./t.^alphat).^(-m),y,1000);
    Er_index = @(y,z)fun_index(y,z);
    Func = @(x,y,z) (-1)^(n+1)*nchoosek(m,n).*y.*x./(1+z)...
        .*exp(-xishu_1(x,y,z)).*exp(-2*pi*lambda.*Er_index(y,z)).*exp(-pi*lambda.*y.^2);
    T1 = @(y,z)arrayfun(@(Y)quadgk(@(x)Func(x,Y,z),0,RL),y);
    T2 = @(z)arrayfun(@(Z)quadgk(@(y)T1(y,Z),0,1000),z);
    
    FuncAsy = @(x,y,z) (-1)^(n+1)*nchoosek(m,n).*y.*x./(1+z)...
        .*(1-xishu_1(x,y,z)).*exp(-2*pi*lambda.*Er_index(y,z)).*exp(-pi*lambda.*y.^2);
    T3 = @(y,z)arrayfun(@(Y)quadgk(@(x)FuncAsy(x,Y,z),0,RL),y);
    T4 = @(z)arrayfun(@(Z)quadgk(@(y)T3(y,Z),0,1000),z);

    index_1(n) = 4*pi*lambda/log(2)/(RL^2)*quadgk(@(z)T2(z),a_low*Gamma,1000);
    index_3(n) = 4*pi*lambda/log(2)/(RL^2)*quadgk(@(z)T4(z),a_low*Gamma,1000);

    
    
    
    xishu_2 = @(x,y) n*miu_t.*Gamma.*Noise.*y.^alphat.*x.^alphat./(P_t*C_RIS_E);
    fun_index2 =@(y) integral(@(t) t-t.*(1+n*miu_t.*y.^alphat./m.*Gamma./t.^alphat).^(-m),y,1000);
    Er_index2 = @(y) arrayfun(@(Y) fun_index2(Y),y);

    index_2(n) =4*pi*lambda/log(2)/(RL^2)*integral3(@(x,y,z) (-1)^(n+1)*nchoosek(m,n).*y.*x./(1+z)...
        .*exp(-xishu_2(x,y)).*exp(-2*pi*lambda.*Er_index2(y)).*exp(-pi*lambda.*y.^2)...
        ,0,RL,0,1000,0,a_low*Gamma);
    index_4(n) =4*pi*lambda/log(2)/(RL^2)*integral3(@(x,y,z) (-1)^(n+1)*nchoosek(m,n).*y.*x./(1+z)...
        .*(1-xishu_2(x,y)).*exp(-2*pi*lambda.*Er_index2(y)).*exp(-pi*lambda.*y.^2)...
        ,0,RL,0,1000,0,a_low*Gamma);
    
    toc;
    
    index_sum(n) = index_1(n) +index_2(n);
    index_sum_Asy(n) = index_3(n) +index_4(n);
    
%     index_sum(n) = index_1(n);
end
CP_t(ii) = sum(index_sum);
CP_t_Asy(ii) = sum(index_sum_Asy);

end

plot(L,CP_t,'r-o');
hold on
plot(L,CP_t_Asy,'b--*');
hold on
