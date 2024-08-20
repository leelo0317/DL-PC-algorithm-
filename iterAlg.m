function [FRF,FBB,WRF,WBB]=iterAlg(H,NRF,Ns)
[U,S,V]=svd(H);
[Nr,Nt]=size(H);
Nbt=Nt/NRF;
Nbr=Nr/NRF;
FRF=[];
WRF=[];
for i = 1:NRF
    FRF = blkdiag(FRF, 1/sqrt(Nbt)*exp(1i * angle(V(Nbt*i-Nbt+1:Nbt*i,i))));
    WRF = blkdiag(WRF, 1/sqrt(Nbr)*exp(1i * angle(U(Nbr*i-Nbr+1:Nbr*i,i))));
end
Heff=WRF'*H*FRF;
Heffinit=zeros(NRF,NRF);
while norm(Heff-Heffinit,'fro')>1e-3
    P=H'*WRF*WRF'*H;
    for m=1:Nt
        c=floor((m-1)/Nbt)+1;
        tep=P([1:m-1, m+1:Nt],m);
        tef=FRF([1:m-1, m+1:Nt],c);
        FRF(m,c)=1/sqrt(Nbt)*exp(1i*angle(tep'*tef));
    end
    Q=H*FRF*FRF'*H';
    for n=1:Nr
        c=floor((n-1)/Nbr)+1;
        teq=Q([1:n-1, n+1:Nr],n);
        tew=WRF([1:n-1, n+1:Nr],c);
        WRF(n,c)=1/sqrt(Nbr)*exp(1i*angle(teq'*tew));
    end
    Heffinit=Heff;
    Heff=WRF'*H*FRF;
    [Ue,Se,Ve]=svd(Heff);
    FBB=Ve(:,1:Ns);
    WBB=Ue(:,1:Ns);
    
end
