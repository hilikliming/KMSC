function [ d_Y ] = KMSC( D, Y, ktype, rbf_var)
% KMSC -- This is a script which implements the M-ary Kernelized Matched
% Subspace Classifier
% Author: John J. Hall (Jack)
%% Inputs:
% D -- Struct containing M matrices of training samples for the M classes
%    - LAST D should always be background samples
%    - e.g. D(1).D = matrix of vector samples from class 1
%
% Y -- Matrix of dimension NxQ observation samples to be classified

%% Outputs:
% d_Y -- matrix of dimension QxM, these are the statistics of each of the Q
% samples under M different hypotheses

d_Y = zeros(size(Y,2),length(D));
for m = length(D)-1
    T_m = D(m).D;           % Class m samples
    Bk  = D(length(D)).D;   % Background samples
    
    KTmTm = kernel(T_m,T_m,ktype,rbf_var);
    KTmB  = kernel(T_m,Bk,ktype,rbf_var);
    KBTm  = kernel(Bk,T_m,ktype,rbf_var);
    KBB   = kernel(Bk,Bk,ktype,rbf_var);
    
    KTBTB = kernel([Tm B],[Tm B],ktype,rbf_var);
    
    [del_m,~,~] = svd(KTBTB);
    
    [T,eT,~] = svd(KTmTm, 'econ');
    [B,eB,~] = svd(KBB, 'econ');
    for nT = 1:size(eT,1)
        if(trace(eT(1:nT,1:nT).^2)/trace(eT.^2)>0.9) 
            break
        end
    end
    for nB = 1:size(eB,1)
        if(trace(eB(1:nB,1:nB).^2)/trace(eB.^2)>0.9) 
            break
        end
    end
    T = T(:,1:nT);
    B = B(:,1:nB);
    
    Lm = [ T'*KTmTm*T    T'*KTmB*B ; ...
           B'*KBTm*T      B'*KBB*B ];
    
    for q = 1:size(Y,2)
        y = Y(:,q);
        % For each sample y we need to calculate:
        KTmBy = kernel([Tm Bk],y, ktype,rbf_var);
        KTmy  = kernel(Tm,y, ktype,rbf_var);
        KBy   = kernel(Bk,y, ktype,rbf_var);
        
        d_Y(q,m) = (norm(del_m'*KTmBy,'fro')^2-norm(B'*KBy,'fro')^2)...
            /(norm(del_m'*KTmBy,'fro')^2-[KTmy'*T KBy'*B]*Lm^(-1)*[T'*KTmy; B'*KBy]);
    end
    
end

end

