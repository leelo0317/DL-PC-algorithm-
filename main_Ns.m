
%%本仿真是大论文第三章（第一篇小论文）对比算法的  数据流与频谱效率和能量效率
%close all
clear
Ncl = 5; % 5; % # of clusters
Nray = 10; % # of rays in each cluster
Nt = 196; %144; % # of transmit antennas
Nr = 36; % # of receive antennas
Nrf = 7; % # of RF chains
Ns = [2,3,4,5,6];
smax = length(Ns);
%realization = 1000;

Pt = 1;
Prf = 0.25;
Pps = 0.05;
Ps = 0.005;

SNR_dB =-5;
SNR = 10.^(SNR_dB./10);
% Pps_2 = 0.014; % # of 2-bit phase shifter
% Pps_3 = 0.015; % # of 3-bit phase shifter
% Pps_4 = 0.045; % # of 3-bit phase shifter

%HNPtest2D196_36NRF7Ns5下的频谱效率记录
% OPT = [16.7861, 22.5300, 26.9713, 30.5644,33.3052];
% SDR = [14.8150, 17.6863, 18.8814, 19.5618, 19.9480 ];
% SIC = [7.2892, 9.6743, 11.2267, 12.5381, 13.4095 ];
%HNPtest2D196_36NRF7Ns4下的频谱效率记录
% OPT = [ ];
% SDR = [  ];
% SIC = [  ];

%注意36_36的比较特殊，用选取的50个训练数据画图，测试数据因类型转换为了单精度不能求解。
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D196_36NRF7Ns4']);
H01 = permute(HNPtest2D196_36NRF7Ns4,[2,3,1]);
H = H01%(:,:,1:50);
%H = HCdb144_36NRFNs4;
realization = size(H,3);

t1 = zeros(smax,realization);
t2 = zeros(smax,realization);

R=zeros(smax,realization);
R1=zeros(smax,realization);
R5=zeros(smax,realization);
RS=zeros(smax,realization);

E=zeros(smax,realization);
E1=zeros(smax,realization);
E5=zeros(smax,realization);
ES=zeros(smax,realization);
for s = 1:smax
s
% for angle_sigma =angle_sigmas
% angle_sigma
% angle_sigma =10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
% gamma = sqrt((Nt*Nr)/(Ncl*Nray)); %normalization factor
% sigma = 1; %according to the normalization condition of the H
% H = zeros(Nr,Nt,realization);
% At = zeros(Nt,Ncl*Nray,realization);
% Ar = zeros(Nr,Ncl*Nray,realization);
% alpha = zeros(Ncl*Nray,realization);

Fopt=zeros(Nt,Ns(s),realization);
Wopt=zeros(Nr,Ns(s),realization);
for reali = 1:realization
%     for c = 1:Ncl
%         AoD_m = unifrnd(0,2*pi,1,2);
%         AoA_m = unifrnd(0,2*pi,1,2);
%         AoD(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(1),angle_sigma);
%         AoD(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(2),angle_sigma);
%         AoA(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(1),angle_sigma);
%         AoA(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(2),angle_sigma);
%     end
%     
%     for j = 1:Ncl*Nray
%         At(:,j,reali) = array_response(AoD(1,j),AoD(2,j),Nt); %UPA array response
%         Ar(:,j,reali) = array_response(AoA(1,j),AoA(2,j),Nr);
%         alpha(j,reali) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
%         H(:,:,reali) = H(:,:,reali) + alpha(j,reali) * Ar(:,j,reali) * At(:,j,reali)';
%     end
%     H(:,:,reali) = gamma * H(:,:,reali);
    if(rank(H(:,:,reali))>=Ns(s))
        [U,S,V] = svd(H(:,:,reali));
        Fopt(:,:,reali) = V([1:Nt],[1:Ns(s)]);
        Wopt(:,:,reali) = U([1:Nr],[1:Ns(s)]);
    end
end


%%
    for reali = 1:realization %par
        reali
%         [FRF1, FBB1] = OMP(Fopt(:,:,reali), Nrf(s), At(:,:,reali));
%         FBB1 = sqrt(Ns(s)) * FBB1 / norm(FRF1 * FBB1,'fro');
%         [WRF1, WBB1] = OMP(Wopt(:,:,reali), Nrf(s), Ar(:,:,reali));       

        %因为SDR算法（部分连接结构）必须要求Nt(Nr)/Nrf是整数，下面这个代码是若Nt(Nr)/Nrf不是整数，则跳过后面的代码，继续进行下次循环
%         if rem(Nt,Nrf(s))~=0 || rem(Nr,Nrf(s))~=0
%             continue
%         end    
        %OPT算法
        R(s,reali) = log2(det(eye(Ns(s)) + SNR/Ns(s) * pinv(Wopt(:,:,reali)) * H(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H(:,:,reali)' * Wopt(:,:,reali)));
        E(s,reali) = R(s,reali)/(Pt+Nt*Prf);
        %SIC算法
        tic
        [ FRFS, FBBS ] = SIC( Fopt(:,:,reali), H(:,:,reali), Nrf, SNR );
         FBBS = sqrt(Ns(s))*FBBS/norm(FRFS*FBBS,'fro');
         
        %toc
        t1(s,reali) = toc;
        
        tic
        %SDR_AltMin算法
        [ FRF_SDR, FBB_SDR ] = SDR_AltMin( Fopt(:,:,reali), Nrf);
        FBB_SDR = sqrt(Ns(s))*FBB_SDR/norm(FRF_SDR*FBB_SDR,'fro');
        %[ WRFS, WBBS ] = SDR_AltMin( Wopt(:,:,reali), Nrf(s));
        t2(s,reali) =  toc;
        %SDR_AltMin算法
        %RS(s,reali)= log2(det(eye(Ns(s)) + SNR/Ns(s) * pinv(WRFS * WBBS) * H(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H(:,:,reali)' * WRFS * WBBS));
        %SDR_AltMin算法 接收端最优
        R_SDR(s,reali)= log2(det(eye(Ns(s)) + SNR/Ns(s) * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRF_SDR * FBB_SDR * FBB_SDR' * FRF_SDR' * H(:,:,reali)' * Wopt(:,:,reali)));
        E_SDR(s,reali) = R_SDR(s,reali)/(Pt+Nrf*Prf+Nt*Pps);
        
        %SIC接收端最优
        RS(s,reali) = log2(det(eye(Ns(s)) + SNR/Ns(s) * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H(:,:,reali)' * Wopt(:,:,reali)));
        ES(s,reali) = RS(s,reali)/(Pt+Nrf*Prf+Nt*Pps);
        %         eS(s,reali)=BER(FRFS,FBBS,WRFS,WBBS,H(:,:,reali),b,SNR_dB);
    end
end

%% 频谱效率
%OPT算法
R_OPT=abs(sum(R,2)/realization);
%FMO算法
%R_FMO=abs(sum(R1,2)/realization);
%HP-HEDC算法
% R2s=abs(sum(R2,2)/realization);
% R3s=abs(sum(R3,2)/realization);
% R4s=abs(sum(R4,2)/realization);
%SIC
RS1=abs(sum(RS,2)/realization);
%SDR_AltMin算法
R_SDR1=abs(sum(R_SDR,2)/realization);
%% 能量效率
E_OPT=abs(sum(E,2)/realization);
%E_FMO=abs(sum(E1,2)/realization);
%HP-HEDC算法
% E2s=abs(sum(E2,2)/realization);
% E3s=abs(sum(E3,2)/realization);
% E4s=abs(sum(E4,2)/realization);

ES1=abs(sum(ES,2)/realization);
E_SDR1=abs(sum(E_SDR,2)/realization);
% %
% es=abs(sum(e,2)/realization);
% e1s=abs(sum(e1,2)/realization);
% e3s=abs(sum(e3,2)/realization);
% e5s=abs(sum(e5,2)/realization);
% eSs=abs(sum(eS,2)/realization);

%对Nt(Nr)/Nrf不是整数的那个点，不画
% Nrf1 = Nrf;
% empty = find(RS1==0);
% Nrf1(empty) = [];
% RS1(empty) = [];
% ES1(empty) = [];
% 
% R_SDR1(empty) = [];
% E_SDR1(empty) = [];
% 
% R_OPT(empty) = [];
% E_OPT(empty) = [];

% Nrf1 = Nrf;
% empty = find(R_SDR1==0);
% Nrf1(empty) = [];
% R_SDR1(empty) = [];
% E_SDR1(empty) = [];
% eSs(s) = [];
%%
figure
plot(Ns, R_OPT,'b--','Marker','+','LineWidth',1.0);
set(gca, 'xtick',2:6);
hold on
grid on
%plot(Nrf, R_FMO,'k--','Marker','s','LineWidth',1.0);

% plot(Nrf, R2s,'y-','Marker','o','LineWidth',1.0);
% plot(Nrf, R3s,'r-','Marker','o','LineWidth',1.0);
% plot(Nrf, R4s,'c-','Marker','o','LineWidth',1.0);

plot(Ns, RS1,'m--','Marker','^','LineWidth',1.0);
plot(Ns, R_SDR1,'g--','Marker','>','LineWidth',1.0);
xlabel('数据流数')%'N_{RF}'
ylabel('频谱效率 (bits/s/Hz)')%'Spectral Efficiency (bits/s/Hz)'
legend('OPT','SIC','SDR-AltMin');%'OMP(FC)','HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)','FMO (FC)',



figure
plot(Ns, E_OPT,'b--','Marker','+','LineWidth',1.0);
set(gca, 'xtick',2:6);
hold on
grid on
%plot(Nrf, E_FMO,'k--','Marker','s','LineWidth',1.0);

% plot(Nrf, E2s,'y-','Marker','o','LineWidth',1.0);
% plot(Nrf, E3s,'r-','Marker','o','LineWidth',1.0);
% plot(Nrf, E4s,'c-','Marker','o','LineWidth',1.0);

plot(Ns, ES1,'m--','Marker','^','LineWidth',1.0);
plot(Ns, E_SDR1,'g--','Marker','>','LineWidth',1.0);
xlabel('数据流数')
ylabel('能量效率 (bits/Hz/J)')%'Energy Efficiency (bits/Hz/J)
legend('OPT','SIC','SDR-AltMin');
%不画全数字
%legend('OMP(FC)','HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)','AS-C-based','SDR-AltMin(PC)');

