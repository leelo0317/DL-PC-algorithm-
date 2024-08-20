
%本仿真是大论文第三章（第一篇小论文）对比算法的 发射天线数与频谱效率和能量效率
clear
%close all
Ns =4; % # of data streams
Ncl = 5; % # of clusters
Nray =10; % # of rays in each cluster
Nr =36; % # of receive antennas
Nrf=4; % # of RF chains

Nt = [36 64 100 144 196]; % # of transmit antennas
smax = length(Nt);

%注意每个H的变量名不同
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D36_36NRF4Ns4']);
H01 = permute(HNPtest2D36_36NRF4Ns4,[2,3,1]);
H1 = H01(:,:,1:50);

load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D64_36NRF4Ns4']);
H02 = permute(HNPtest2D64_36NRF4Ns4,[2,3,1]);
H2 = H02(:,:,1:50);
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D100_36NRF4Ns4']);
H03 = permute(HNPtest2D100_36NRF4Ns4,[2,3,1]);
H3 = H03(:,:,1:50);
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D144_36NRF4Ns4']);
H04 = permute(HNPtest2D144_36NRF4Ns4,[2,3,1]);
H4 = H04(:,:,1:50);
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D196_36NRF4Ns4']);
H05 = permute(HNPtest2D196_36NRF4Ns4,[2,3,1]);
H5 = H05(:,:,1:50);
% 将这些信道存储在H元胞数组中，以便下面循环访问
Hcell = {H1, H2, H3, H4, H5}
%下面循环访问中，求Fopt时会用到cell的访问，规则如下： 参看网址：
% cell的访问:
% 1.X= C(s)使用这种"()"形式的返回的是cell类
% 2.X = C{s}使用这种"{}"形式的返回的是cell中的内容
% 参看网址：https://blog.csdn.net/u011171235/article/details/51148519

%50为信道矩阵的个数
Fopt1 = zeros(36,4,100);
Fopt2 = zeros(64,4,100);
Fopt3 = zeros(100,4,100);
Fopt4 = zeros(144,4,100);
Fopt5 = zeros(196,4,100);
Foptcell = { Fopt1, Fopt2, Fopt3, Fopt4, Fopt5};
realization = 100;
% i=1;
Pt = 1;
Prf = 0.25;
Pps = 0.05;
Ps = 0.005;
%原刘华菁HEDC代码
% Pps_2 = 0.014; % # of 2-bit phase shifter
% Pps_3 = 0.015; % # of 3-bit phase shifter
% Pps_4 = 0.045; % # of 3-bit phase shifter
%为频谱效率预先分配空间
R=zeros(smax,realization);
%R_FHB=zeros(smax,realization);
% R2=zeros(smax,realization);
% R3=zeros(smax,realization);
% R4=zeros(smax,realization);
RSIC=zeros(smax,realization);
RS=zeros(smax,realization);
REGT=zeros(smax,realization);
%为能量效率预先分配空间
E=zeros(smax,realization);
%E_FHB=zeros(smax,realization);
%
% E2=zeros(smax,realization);
% E3=zeros(smax,realization);
% E4=zeros(smax,realization);
ESIC=zeros(smax,realization);
ES=zeros(smax,realization);
EEGT=zeros(smax,realization);
% e=zeros(smax,realization);
% e1=zeros(smax,realization);
% e3=zeros(smax,realization);
% e5=zeros(smax,realization);
% eS=zeros(smax,realization);

for s = 1:smax
    s
%     angle_sigma =10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
%     gamma = sqrt((Nt(s)*Nr)/(Ncl*Nray)); %normalization factor
%     sigma = 1; %according to the normalization condition of the H
%     H = zeros(Nr,Nt(s),realization);
%     At = zeros(Nt(s),Ncl*Nray,realization);
%     Ar = zeros(Nr,Ncl*Nray,realization);
%     alpha = zeros(Ncl*Nray,realization);
%     Fopt=zeros(Nt(s),Ns,realization);
%     Wopt=zeros(Nr,Ns,realization);
%     for reali = 1:realization
%         for c = 1:Ncl
%             AoD_m = unifrnd(0,2*pi,1,2);
%             AoA_m = unifrnd(0,2*pi,1,2);
%             AoD(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(1),angle_sigma);
%             AoD(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(2),angle_sigma);
%             AoA(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(1),angle_sigma);
%             AoA(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(2),angle_sigma);
%         end
%     
%         for j = 1:Ncl*Nray
%             At(:,j,reali) = array_response(AoD(1,j),AoD(2,j),Nt(s)); %UPA array response
%             Ar(:,j,reali) = array_response(AoA(1,j),AoA(2,j),Nr);
%             alpha(j,reali) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
%             H(:,:,reali) = H(:,:,reali) + alpha(j,reali) * Ar(:,j,reali) * At(:,j,reali)';
%         end
%         H(:,:,reali) = gamma * H(:,:,reali);
    
    H = Hcell{s};
    Fopt = Foptcell{s};
    
    realization = size(H,3);
    for reali = 1:realization
        if(rank(H(:,:,reali))>=Ns)
            [U,S,V] = svd(H(:,:,reali));
            Fopt(:,:,reali) = V([1:Nt(s)],[1:Ns]);
            Wopt(:,:,reali) = U([1:Nr],[1:Ns]);
        end
    end

    SNR_dB =15;
    SNR = 10.^(SNR_dB./10);

    %b=8;
%%
    parfor reali = 1:realization
        reali
        %
%         [FRF1, FBB1] = OMP(Fopt(:,:,reali), Nrf, At(:,:,reali));
%         FBB1 = sqrt(Ns) * FBB1 / norm(FRF1 * FBB1,'fro');
%         [WRF1, WBB1] = OMP(Wopt(:,:,reali), Nrf, Ar(:,:,reali));
        %     %
%         %PE_AltMin算法
%         [ FRF, FBB ] = PE_AltMin( Fopt(:,:,reali), Nrf);
%         FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
%         [ WRF, WBB ] = PE_AltMin( Wopt(:,:,reali), Nrf);   
        tic
        %SDR_AltMin算法
        [ FRFS, FBBS ] = SDR_AltMin( Fopt(:,:,reali), Nrf);
        FBBS = sqrt(Ns)*FBBS/norm(FRFS*FBBS,'fro');
        %[ WRFS, WBBS ] = SDR_AltMin( Wopt(:,:,reali), Nrf);
        toc      
        %SIC
        tic
        [ FRFSIC, FBBSIC ] = SIC( Fopt(:,:,reali), H(:,:,reali), Nrf, SNR );
        FBBSIC = sqrt(Ns)*FBBSIC/norm(FRFSIC*FBBSIC,'fro');
        toc
       
        tic
        [FRF7,FBB7,WRF7,WBB7]=iterAlg(H(:,:,reali),Nrf,Ns);%等增益（白晓宇大论文第3章）
        FBB7 = sqrt(Ns)*FBB7/norm(FRF7*FBB7,'fro');
        toc
        
        %EGT
        REGT(s,reali)= log2(det(eye(Ns) + SNR/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRF7 *FBB7* (FRF7*FBB7)' * H(:,:,reali)' * Wopt(:,:,reali)));
        EEGT(s,reali) = REGT(s,reali)/(Pt+Nrf*Prf+Nt(s)*Pps);
        
        %最优OPT
        R(s,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H(:,:,reali)' * Wopt(:,:,reali)));
        E(s,reali) = R(s,reali)/(Pt+Nt(s)*Prf);
%         e(s,reali)=BER(Fopt(:,:,reali),eye(Ns),Wopt(:,:,reali),eye(Ns),H(:,:,reali),b,SNR_dB);
        %OMP算法
%         R1(s,reali)= log2(det(eye(Ns) + SNR/Ns * pinv(WRF1 * WBB1) * H(:,:,reali) * FRF1 * FBB1 * FBB1' * FRF1' * H(:,:,reali)' * WRF1 * WBB1));
%         E1(s,reali) = R1(s,reali)/(Pt+Nrf*Prf+Nrf*Nt(s)*Pps);
        %FHB算法
%         R_FHB(s,reali)= log2(det(eye(Ns) + SNR/Ns * pinv(WRFk * WBBk) * H(:,:,reali) * FRFk * FBBk * FBBk' * FRFk' * H(:,:,reali)' * WRFk * WBBk));
%         E_FHB(s,reali) = R_FHB(s,reali)/(Pt+Nrf*Prf+Nrf*Nt(s)*Pps);      
        %SIC
        RSIC(s,reali)= log2(det(eye(Ns) + SNR/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRFSIC *FBBSIC * (FRFSIC*FBBSIC)' * H(:,:,reali)' * Wopt(:,:,reali)));
        ESIC(s,reali) = RSIC(s,reali)/(Pt+Nrf*Prf+Nt(s)*Pps);
%         e5(s,reali)=BER(FRF5 *FS5, FBB5,WRF5 *WS5, WBB5,H(:,:,reali),b,SNR_dB);
        
        %SDR_AltMin算法
        RS(s,reali)= log2(det(eye(Ns) + SNR/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H(:,:,reali)' * Wopt(:,:,reali)));
        ES(s,reali) = RS(s,reali)/(Pt+Nrf*Prf+Nt(s)*Pps);
%         eS(s,reali)=BER(FRFS,FBBS,WRFS,WBBS,H(:,:,reali),b,SNR_dB);

    end
end
R_opt=abs(sum(R,2)/realization);
%R_FHBs=abs(sum(R_FHB,2)/realization);
R_SIC=abs(sum(RSIC,2)/realization);
R_SDR=abs(sum(RS,2)/realization);
R_EGT=abs(sum(REGT,2)/realization);
%
E_opt=abs(sum(E,2)/realization);
%E_FHBs=abs(sum(E_FHB,2)/realization);
E_SIC=abs(sum(ESIC,2)/realization);
E_SDR=abs(sum(ES,2)/realization);
E_EGT=abs(sum(EEGT,2)/realization);
% %
% es=abs(sum(e,2)/realization);
% e1s=abs(sum(e1,2)/realization);
% e3s=abs(sum(e3,2)/realization);
% e5s=abs(sum(e5,2)/realization);
% eSs=abs(sum(eS,2)/realization);
figure()
%plot(Nt, Rs,'b--','Marker','+','LineWidth',1.0);
axis([36 196 20 55]);
set(gca, 'xtick',[36, 64,100,144,196]);
hold on
grid on
plot(Nt, R_opt,'k--','Marker','s','LineWidth',1.0);
%所提
% plot(Nt, R2s,'y-','Marker','o','LineWidth',1.0);
% plot(Nt, R3s,'r-','Marker','o','LineWidth',1.0);
% plot(Nt, R4s,'c-','Marker','o','LineWidth',1.0);

plot(Nt, R_SIC,'m--','Marker','^','LineWidth',1.0);
plot(Nt, R_SDR,'g--','Marker','>','LineWidth',1.0);
plot(Nt, R_EGT,'c--','Marker','>','LineWidth',1.0);
xlabel('N_t')
ylabel('Spectral Efficiency (bits/s/Hz)')
%legend('Fully Digitial','OMP(FC)','HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)','AC-S-based','SDR-AltMin(PC)');
%不画OPT
legend('OPT','SIC','SDR-AltMin','EGT');%'HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)',




figure()
%plot(Nt, Es,'b--','Marker','+','LineWidth',1.0);
axis([36 196 0 10]);
set(gca, 'xtick',[36, 64,100,144,196]);
hold on
grid on
plot(Nt, E_opt,'k--','Marker','s','LineWidth',1.0);
%所提
% plot(Nt, E2s,'y-','Marker','o','LineWidth',1.0);
% plot(Nt, E3s,'r-','Marker','o','LineWidth',1.0);
% plot(Nt, E4s,'c-','Marker','o','LineWidth',1.0);

plot(Nt, E_SIC,'m--','Marker','^','LineWidth',1.0);
plot(Nt, E_SDR,'g--','Marker','>','LineWidth',1.0);
plot(Nt, E_EGT,'c--','Marker','>','LineWidth',1.0);
xlabel('N_t')
ylabel('Energy Efficiency (bits/Hz/J)')
%legend('Fully Digitial','OMP(FC)','HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)','AC-S-based','SDR-AltMin(PC)');
%不画OPT
legend('OPT','SIC','SDR-AltMin','EGT');%,'HP-HEDC(2bit)','HP-HEDC(3bit)','HP-HEDC(4bit)'
