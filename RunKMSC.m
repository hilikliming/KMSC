%% Kernelized Matched Subspace Classifier
% John Hall
% Based off the work of H. Kwon and N. Nasrabadi in their 'Hyperpsectral 
% target detection using Matched Subspace Detectors in the original input 
% domain'

% L_m(y) = 
%   num: ( K(Z_{T_m,B},y)'*del_m*del_m'*K(Z_{T_m,B},y)-K(Z_{B},y)'*B*B'*K(Z_{B},y)  )./
%   den: ( K(Z_{T_m,B},y)'*del_m*del_m'*K(Z_{T_m,B},y)-[K(Z_{T_m},y)'*T, ...
% K(Z_{B},y)'*B]*Lm^(-1)*[T'*K(Z_{Tm},y ; B'*K(Z_{B},y] )

%this can be further simplfied to:
% L_m(y) = 
% num: ( |del_m'*K(Z_{T_m,B},y)|^2_F -|B'*K(Z_{B},y)|_F^2 )./
% den: ( |K(Z_{T_m,B},y)'*del_m|^2_F - [K(Z_{T_m},y)'*T, ...
%        K(Z_{B},y)'*B]*Lm^(-1)*[T'*K(Z_{Tm},y ; B'*K(Z_{B},y] )


% Here
% Lm = [ T'*K(Z_{T_m},Z_{T_m})*T    T'*K(Z_{T_m},Z_{B})*B ; ...
%        B'*K(Z_{B},Z_{T_m})*T      B'*K(Z_{B},Z_{B})*B ]

%del_m are ALL eigenvectors of K(Z_{T_m,B},Z_{T_m,B}) stacked side by side

% Matrices B and T are the SIGNIFICANT eigenvectors of background centered 
% K(Z_{B},Z_{B}) and K(Z_{T_m},Z_{T_m}) respectively...

% Kc = (K-o_m*K-K*o_m + o_m*K*o_m
% where o_m = ones(size(K))./size(K,1);

rbf_var = 3;
ktype = 'rbf';

% For more compact notation use Z_{T_m} = T_m and Z_{B} = B
% These are not the 'T' and B from above, these are the training samples of
% the T_m and B classes.

% For any sample, to form the likelihood ratio, under mth hypothesis we need:
% 1) K(T_m,T_m)
% 2) K(T_m,  B)
% 3) K(B ,T_m)
% 4) K(B,B)
% 5) [del_m ~, ~] = svd(K([T_m,B],[T_m,B]))
% 6) Lm = [ T'*K(Z_{T_m},Z_{T_m})*T    T'*K(Z_{T_m},Z_{B})*B ; ...
%           B'*K(Z_{B},Z_{T_m})*T      B'*K(Z_{B},Z_{B})*B ]

% For each sample y we need to calculate:
% 1) K([T_m,B],y)
% 2) K(T_m,y)
% 3) K(B,y)


%% Importing Training Data


%% Running KMSC

%% Documenting Results


