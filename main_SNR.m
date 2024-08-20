
%此代码为部分连接对比算法（SDR和SIC算法）
clear,clc
%addpath(pwd);
%cd cvx;
%addpath(genpath(pwd));
%cd ..;

% load('Ns=3');
% H_Train = H(:,:,1:100)

%load(['F:\tiangui_bufen\dataset\5-10-15dB\Ns=NRF=4\H_Train100Ns4']);
load(['F:\tiangui_bufen\dataset\H5_10_15dB\test_dataWopt\HNPtest2D196_36NRF7Ns4']);
H01 = permute(HNPtest2D196_36NRF7Ns4,[2,3,1]);
H_Train = H01;
%load('F:\tiangui_bufen\dataset\perfect_H\imperfect_H_2D');

Nt=196;
Nr=36;
NRF = 7;
Ns = 7;
%H_Train = awgn(H_Train ,-20);
%M=Nt/NRF; % number of antennas connected to one RF chains
SNR_dB = -20:5:20; %-35:5:5;
%SNR_dB = -40:2:-20
SNR = 10.^(SNR_dB./10);
smax = length(SNR);% enable the parallel
realization = size(H_Train,3);
%存储每个信道矩阵带入相应算法中生成混合预编码的时间
t1 = zeros(1,realization);
t2 = zeros(smax,realization);

%%
%若load的信道数据里没有 Fopt和 Wopt，用下面一段代码来生成，以便SDR_AltMins使用
for reali = 1:realization    
    
    if(rank(H_Train(:,:,reali))>=Ns)
        [U,S,V] = svd(H_Train(:,:,reali));
        Fopt(:,:,reali) = V([1:Nt],[1:Ns]);
        
        Wopt(:,:,reali) = U([1:Nr],[1:Ns]);
    end
end
%%

%SNR_dB = -35:5:5;
%SNR = 10.^(SNR_dB./10);
%realization = size(H,3);

%下面的SDR和SIC算法都是假设的接收端最优
for reali = 1:realization
    
    
    reali
    tic 
    [ FRF, FBB ] = SDR_AltMin( Fopt(:,:,reali), NRF);
    FBB = sqrt(Ns)*FBB/norm(FRF*FBB,'fro');
    F1 =  FRF*FBB;
    t1(1,reali)  = toc;
    %[ WRF, WBB ] = Receiver( Wopt(:,:,reali), NRF);
    
    for s = 1:smax

        tic
        [ FRFS, FBBS ] = SIC( Fopt(:,:,reali), H_Train(:,:,reali), NRF, SNR(s) );
        FBBS = sqrt(Ns)*FBBS/norm(FRFS*FBBS,'fro');
        F2 =  FRFS*FBBS;
        t2(s,reali) = toc;
        
        %SIC接收端最优
        RS(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H_Train(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H_Train(:,:,reali)' * Wopt(:,:,reali)));
        
        %SDR
        %R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H_Train(:,:,reali) * FRF * FBB * FBB' * FRF' * H_Train(:,:,reali)' * WRF * WBB));
        %SDR接收端最优
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H_Train(:,:,reali) * FRF * FBB * FBB' * FRF' * H_Train(:,:,reali)' * Wopt(:,:,reali)));
  
        %[ WRFS, WBBS ] = SIC( Wopt(:,:,reali), H_Train(:,:,reali), NRF, SNR(s) );
        
        %改SIC成接收端不是最优,接收端利用SDR
        %RS(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H_Train(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H_Train(:,:,reali)' *WRF * WBB));
        %不考虑接收端时的最优
        %Ro(s,reali) = log2(det(eye(Nr) + SNR(s)/Ns * H_Train(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H_Train(:,:,reali)' ));
  
        %OPT
        Ro(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H_Train(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H_Train(:,:,reali)' * Wopt(:,:,reali)));
        
    end
end
figure
grid on
hold on
plot(SNR_dB,abs(sum(R,2)/realization),'Marker','diamond','LineWidth',1.5,'Color',[0.87058824300766 0.490196079015732 0]);
plot(SNR_dB,abs(sum(RS,2)/realization),'c-+','LineWidth',1.5);
plot(SNR_dB,abs(sum(Ro,2)/realization),'r-o','LineWidth',1.5);
%plot(SNR_dB,abs(sum(R3,2)/realization),'y-s','LineWidth',1.5);
%
xlabel('SNR [dB]');
ylabel('Spectral Efficiency [bits/s/Hz]');
legend('SDR-AltMin','SIC','OPT');%,'Analog'

% A=[5.689381103515625, 9.622802734375, 14.140501708984376, 18.967021484375, 23.91871337890625, 28.91634765625, 33.93334228515625, 38.96603515625, 44.0178662109375];
% % %[1.9469020080566406, 4.05562744140625, 7.19663330078125, 11.2476904296875, 15.96416748046875, 21.13756103515625, 26.6588525390625, 32.4971923828125, 38.631396484375];
% % % %[3.1411776733398438, 6.672532348632813, 11.745726318359376, 17.81698974609375, 24.3436376953125, 31.0393603515625, 37.79190185546875, 44.5627880859375, 51.339541015625]
% plot(SNR_dB,A,'y-s','LineWidth',1.5);
