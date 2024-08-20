function [FRF,FBB] = SDR_AltMin(Fopt,NRF)

% randomly generate FRF
[Nt,Ns] = size(Fopt);
FRF = [];
for i = 1:NRF
    FRF = blkdiag(FRF, exp(sqrt(-1) * unifrnd (0,2*pi,[Nt/NRF,1])));
end
FRF = 1/sqrt(Nt)*FRF;

y = [];
while(isempty(y) || abs(y(1)-y(2))>1e-3)
    % fix FRF, optimize FBB
    A1 = diag([ones(1,Ns*NRF),0]);
    A2(Ns*NRF+1,Ns*NRF+1) = 1;
    
    temp = kron(eye(Ns),FRF);
    C = [temp'*temp,-temp'*vec(Fopt);-vec(Fopt)'*temp,vec(Fopt)'*vec(Fopt)];
    
    cvx_begin quiet
        variable Y(Ns*NRF+1,Ns*NRF+1) hermitian %我把X改成了Y
        minimize(real(trace(C*Y)));%我把X改成了Y
        subject to
            trace(A1*Y) == NRF*Ns;%我把X改成了Y
            trace(A2*Y) == 1;%我把X改成了Y
        Y == hermitian_semidefinite(Ns*NRF+1);%我把X改成了Y
    cvx_end
    
    [V,D] = eig(Y);%我把X改成了Y
    [value,num] = max(diag(D));
    x = sqrt(value)*V(:,num);
    FBB = reshape(x(1:Ns*NRF),NRF,Ns);
    
    y(1) = norm(Fopt-FRF*FBB,'fro')^2;
    
    % fix FBB, optimize FRF
    for i = 1:Nt
        m = ceil(i*NRF/Nt);
        FRF(i,m) = 1/sqrt(Nt) * exp( sqrt(-1) * angle(Fopt(i,:)*FBB(m,:)') );
    end
    
    y(2) = norm(Fopt-FRF*FBB,'fro')^2;
end

end