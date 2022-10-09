%% SETTINGS %%
clear
clc
L = 8; % half length of RIS
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
r_C = 70; % connedted user distance to nearest user

Pt_dBm = -5:5:10; % transmit power dBm
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
%% SIMULATION %%
Pc_cRIS = zeros(1,length(SNR_dB));
Pc_tRIS = zeros(1,length(SNR_dB));
for ii=1:length(SNR_dB)
    tic
SNR_RIS = zeros(1,num);
SNR_c_RIS = zeros(1,num);
SNR_SIC = zeros(1,num);
channel_t = gamrnd(m,1/m,[1,num]); % gamma dist. for typical user
channel_c = gamrnd(m,1/m,[1,num]); % gamma dist. for connec. user
for k=1:num
    
%% distance 
% settings_BR
N_bs=poissrnd(lambda.*pi.*R_BS.^2); % number of BSs
rc=sqrt(R_BS.^2.*rand(N_bs,1)); % BSs in polar coordinates (the distances from a typical user as the center)
trans=find(rc==min(rc)); % the closest UAV is selected as the typical UAV
r_BR=rc(trans); % the distance from nearest BS to RIS
rc(trans)=[]; % the rest distances, these BSs as interferences
r_BR_I = rc; % distances from interference BSs to typical user
%RU
r_RU=sqrt(RL.^2.*rand(1,1)); % RISs in LOS ball



%% angles of incidence and refelection for nearst user
% psi_1 = 2.*pi.*rand(N_bs,1); % user angle from PPP
% psi_2 = 2.*pi.*rand(1,1); % RIS angle from LOS ball
% psi_1_near = psi_1(trans);
% psi_1(trans)=[];
% psi_1_I = psi_1;
% angles of incidence and refelection for nearest BS
theta_psi = pi.*rand(1,1);
% theta_psi = pi.*0.6;
% rho_r = rand(1,1);
rho_r = 0.4;
theta_inc = rho_r*theta_psi; % angles of incidence = abgles of reflection
theta_ref = (1-rho_r)*theta_inc; % abgles of reflection
% angles of incidence and refelection for Interference
theta_psi_I =  pi.*rand(N_bs-1,1);
% theta_psi_I =  pi.*0.6;
rho_r_I =  0.4;
theta_inc_I = rho_r_I.*theta_psi_I; % angles of incidence = abgles of reflection
theta_ref_I = (1-rho_r_I).*theta_inc_I; % abgles of reflection


% %% calculate the intercrept of interference
% rho_t_theta = 2*pi.*rand(length(r_BR_I),1); % the angles between typical user and BSs are uniformly distributed
% rho_t = zeros(length(r_BR_I),1); % store space for rho_t
% for aa = 1:length(r_BR_I)
%     if rho_t_theta(aa)>=pi
%         rho_t(aa) = 0;
%     elseif  rho_t_theta(aa)<pi
%         rho_t(aa) = 1;
%     end
% end

%% Path loss 
% for nearest user
C_RIS_near = L/4/pi.*abs(cos(theta_inc)+cos(theta_ref));
P_t_RIS_near(k) = C_RIS_near^2*(r_BR*r_RU)^(-alphat); % RIS path loss signal
% for interference
C_RIS_I = L/4/pi.*abs(cos(theta_inc_I)+cos(theta_ref_I));
P_t_RIS_I = C_RIS_I.^2.*(r_BR_I.*r_RU).^(-alphat); % RIS path loss interference
% I_t_RIS(k) = sum(rho_t.*P_t(ii).*gamrnd(m,1/m,[length(r_BR_I),1]).*P_t_RIS_I); % total interference of typical user
I_t_RIS(k) = 0.5*sum(P_t(ii).*gamrnd(m,1/m,[length(r_BR_I),1]).*P_t_RIS_I);

%% connected user
r_C_I = sqrt((R_BS.^2-r_C^2).*rand(N_bs-1,1)); % distence for the connected user
P_c = P_t(ii).*C_L*r_C^(-alpha);% connected user path loss
P_c_I = P_t(ii).*C_L.*r_C_I.^(-alpha);% connected user path loss interference 
RP_C_I(k) = sum( gamrnd(m,1/m,[length(P_c_I),1]).* P_c_I) ;   

%% SNR 
% connected
SNR_c_RIS(k) = a_high*channel_c(k)* P_c./ ...
              ( a_low.*channel_c(k)*P_c + RP_C_I(k) + Noise);% Connected User
SNR_c_OMA(k) = channel_c(k)* P_c./( RP_C_I(k) + Noise);% OMA
% typical
SNR_SIC(k) = a_high*P_t(ii)*channel_t(k)* P_t_RIS_near(k) ./...
             ( a_low.*P_t(ii)*channel_t(k)* P_t_RIS_near(k) + I_t_RIS(k) + Noise);
SNR_RIS(k)= a_low.*P_t(ii)*channel_t(k)* P_t_RIS_near(k)/ ( I_t_RIS(k) + Noise); % SINR Typical User 
SNR_t_OMA(k) = P_t(ii)*channel_t(k)* P_t_RIS_near(k)./( I_t_RIS(k) + Noise); %OMA
SNR_t_OMA_Noise0(k) = P_t(ii)*channel_t(k)* P_t_RIS_near(k)./( I_t_RIS(k) ); %OMA
end

% % OMA
% threshold_c(ii) =  mean(channel_c)* P_c./(mean(RP_C_I) + Noise);
threshold_c2(ii) = mean(channel_c.* P_c./(RP_C_I + Noise));
% threshold_t(ii) = P_t(ii)*mean(channel_t)* mean(P_t_RIS_near)./(mean(I_t_RIS) + Noise);
threshold_t2(ii) = mean(P_t(ii).*channel_t.* P_t_RIS_near./(I_t_RIS + Noise));

%% Ergodic Rate

% Pc_cRIS(ii) = sum( SNR_c_OMA<threshold_t2(ii) )/num;
% Pc_tRIS(ii)=sum( SNR_SIC>threshold2  & threshold_c2(ii)< SNR_t_OMA )/num;
% Er_c(ii) = mean(log2(1+SNR_c_RIS))*Pc_cRIS(ii);
% Er_t(ii) = mean(log2(1+SNR_RIS))*Pc_tRIS(ii);
for  jj = 1:num
    if SNR_c_OMA(jj)<threshold_t2(ii)
        Er_c_index(jj)= mean(log2(1+SNR_c_RIS));
    else 
        Er_c_index(jj)=0;
    end 
    if SNR_SIC(jj)>threshold2  && threshold_c2(ii)< SNR_t_OMA(jj)
        Er_t_index(jj) = mean(log2(1+SNR_RIS));
    else 
        Er_t_index(jj) = 0;
    end
end
   Er_c(ii) = mean(Er_c_index)     ;
   Er_t(ii) = mean(Er_t_index);
toc
end

%% FIGHURE %%
% % figure
% plot(SNR_dB,Er_c,'k-*');
% hold on
plot(SNR_dB,Er_t,'k-o');
hold on
xlabel('Transmit SNR $\rho=P_t/\sigma^2$ (dB)');
ylabel('Ergodic Rate (BPCU)');
legend('1','2');


