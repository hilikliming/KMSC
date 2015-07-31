function [ d_Y ] = KMSC( D, Y, ktype, rbf_var,per)
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

d_Y = zeros(size(Y,2),length(D)-1);
Ztot=[];

% Collecting all training samples 
for i = 1:length(D)
    Ztot = [Ztot , D(i).D];
end

KTBTB               = kernel(Ztot,Ztot,ktype,rbf_var);
[del_tot,deltot,~]  = svd(KTBTB);
del_tot = del_tot*(deltot^(-1/2)); %dividing out square root of eigenvalue

Bk  = D(length(D)).D;   % Background samples

for m = 1:length(D)-1
    T_m = D(m).D;           % Class m samples
    
    % Creating Kernel Matrices for Target Class M and Background Samples
    KTmTm = kernel(T_m,T_m,ktype,rbf_var);
    KTmB  = kernel(T_m,Bk,ktype,rbf_var);
    KBTm  = kernel(Bk,T_m,ktype,rbf_var);
    KBB   = kernel(Bk,Bk,ktype,rbf_var);
    
    % Taking Eigenvectors from KTmTm and KBB
    [T,eT,~] = svd(KTmTm, 'econ');
    [B,eB,~] = svd(KBB, 'econ');
    
    %nBT =size(del_tot,2);
%         if(trace(delm(1:nBT,1:nBT))/trace(delm)>=1) 
%             break
%         end
%     end
    for nT = 1:size(eT,2)
        if(trace(eT(1:nT,1:nT))/trace(eT)>per) 
            break
        end
    end
    for nB = 1:size(eB,2)
        if(trace(eB(1:nB,1:nB))/trace(eB)>per) 
            break
        end
    end
    %Only taking top Eigenvectors
    
    T       = T*(eT^(-1/2));
    B       = B*(eB^(-1/2));
    
    %del_tot = del_tot(:,1:nBT);
    T = T(:,1:nT);
    B = B(:,1:nB);
    
    Lm = [ T'*KTmTm*T    T'*KTmB*B ; ...
           B'*KBTm*T      B'*KBB*B ];
    
    KTBY    = kernel(Ztot,Y, ktype,rbf_var); 
    KTmY    = kernel(T_m,Y, ktype,rbf_var); 
    KBY     = kernel(Bk,Y, ktype,rbf_var); 
    
    for q = 1:size(Y,2)
        % For each sample y we need to calculate:
        KTBy  = KTBY(:,q);
        KTmy  = KTmY(:,q);
        KBy   = KBY(:,q);
        
         %(norm(del_m'*KTmBy,'fro')^2-norm(B'*KBy,'fro')^2);%...
        d_Y(q,m) =(norm(del_tot'*KTBy,'fro')^2-[KTmy'*T KBy'*B]*(Lm^(-1))*[T'*KTmy; B'*KBy]);%/kernel(Y(:,q),Y(:,q),ktype,rbf_var);%/trace(KTBy'*KTBy);
        %warning('off','last');
    end
    
end
d_Y=normc(d_Y);
end

