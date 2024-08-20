function [ FRF, FBB ] = OMP(Fopt, NRF, At)
Ns=size(Fopt,2);
FRF = [];
Fres = Fopt;
for k = 1:NRF
    PU = At' * Fres;
%     [aa,bb] = max(diag(PU * PU'));
    [aa,bb] = sort(sum( abs(PU).^2, 2 ),'descend');
    FRF = [FRF , At(:,bb(1))];
    FBB = pinv(FRF) * Fopt; %use pseudoinverse to avoid the inverse of a possible singular matrix
    Fres = (Fopt - FRF * FBB) / norm(Fopt - FRF * FBB,'fro');
end

end

