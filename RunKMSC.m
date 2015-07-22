home = cd;
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

% del_m are ALL eigenvectors of K(Z_{T_m,B},Z_{T_m,B}) stacked side by side

% Matrices B and T are the SIGNIFICANT eigenvectors of background centered 
% K(Z_{B},Z_{B}) and K(Z_{T_m},Z_{T_m}) respectively...

% Kc = (K-o_m*K-K*o_m + o_m*K*o_m
% where o_m = ones(size(K))./size(K,1);

% For more compact notation use Z_{T_m} = T_m and Z_{B} = B
% These are not the 'T' and B from above, these are the training samples of
% the T_m and B classes.

% For any sample, to form the likelihood ratio, under mth hypothesis we need:
% 1) K(T_m,T_m)
% 2) K(T_m,  B)
% 3) K(B  ,T_m)
% 4) K(B  ,  B)
% 5) [del_m ~, ~] = svd(K([T_m,B],[T_m,B]))
% 6) Lm = [ T'*K(Z_{T_m},Z_{T_m})*T    T'*K(Z_{T_m},Z_{B})*B ; ...
%           B'*K(Z_{B},Z_{T_m})*T      B'*K(Z_{B},Z_{B})*B ]

% For each sample y we need to calculate:
% 1) K([T_m,B],y)
% 2) K(T_m,y)
% 3) K(B,y)


%% Importing Training Data
dirm = 'C:\Users\halljj2\Desktop\TESTMATS';
above = '..\';
cd(above);
dirFRM  = 'DBFRM';
cd(dirFRM); 
dirFRM = cd;

% % Parameters for generateDatabaseLsas indicates environment conditions,
% % ranges, and target rotations to be modeled
% ranges = 10:5:40;
% %water and sediment sound speeds
% c_w = [1464,1530];
% c_s = [1694,1694];
% % rotations to model
% rots = [0:20:80,270:20:350];
% % environment parameters to model, water, sediment speed, interface elevation
% envs      = zeros(length(c_w),3);
% envs(:,1) = c_w;
% envs(:,2) = c_s;
% envs(:,3) = 3.8*ones(length(c_w),1);
% % which of the 7 .ffn's to model
% objs    = 1:10; 
% f_s     = 1e5;
% 
% chirp.sigName  = 'TREXred.replica';%
% chirp.chirp = [1 31 f_s]; % start and end freq of chirp defines center and BW, last number is f_s
% runlen  = [20,800]; %length in meters, stops
% 
% cd(home);
% dirMapDBFRM = generateDatabaseLsas(dirFRM,envs,ranges,rots,objs,chirp,runlen);
% dirm = 'C:\Users\halljj2\Desktop\TESTMATS';
 cd(dirm);
% save('dirMapDBFRM.mat','dirMapDBFRM');

open('dirMapDBFRM.mat');
cd(home);

%% Collecting TRAINING Data
% NT = [1 2 3];
% T = [4 5 6 7 8 9 10];
% stops = 800;
% D = struct([]);
% 
% D(1).D = getTrainingSamples(dirMapDBFRM(:,NT,:,:), stops);
% D(2).D = getTrainingSamples(dirMapDBFRM(:, T,:,:), stops);
% %D(3).D = getTrainingSamples(dirMapDBFRM(:,BK,:,:), stops);
% cd(dirm);
% save('D.mat','D')
cd(dirm);
open('D.mat');
cd(home);
%% Importing Testing Data
% realTarg = {'TARGET17','TARGET7','TARGET16','TARGET20','TARGET21',...
%     'TARGET25','TARGET29','TARGET9','TARGET28',...
%     'TARGET8'}; % 17,7,16 are NT, rest are T
% 
% rngs = [0,10,0,0,0,0,0,0;    ...% 1)  Target 17 alcyl_2ft      NT
%         0,0,0,0,0,30,35,40;  ...% 2)  Target 7  alcyl_3ft      NT
%         0,0,15,0,25,30,0,0;  ...% 3)  Target 16 alpipe         NT
%         0,10,15,0,0,30,0,40; ...% 4)  Target 20 aluxo          T
%         0,10,15,0,25,30,0,0; ...% 5)  Target 21 ssuxo          T
%         5,0,15,20,25,0,0,0; ... % 6)  Target 25 bullet105air   T
%         0,0,15,20,0,0,35,0;  ...% 7)  Target 29 bullet105h20   T
%         0,10,15,20,0,0,35,0; ...% 8)  Target 9  howcapair      T
%         0,0,15,0,25,30,0,0;  ...% 9)  Target 28 howcaph20      T
%         5,0,0,0,25,30,0,40];    % 10) Target 8  hownocap        T   
% eps = 110;
% [Y, t_Y, r_Y ] = realACfetchTREX(realTarg,stops,eps,rngs);
% cd(dirm)
% save('Y.mat','Y');
% save('t_Y.mat','t_Y');
% save('r_Y.mat','r_Y');

cd(dirm);
open('Y.mat');
open('r_Y.mat');
open('t_Y.mat');
cd(home);
%% Running KMSC
rbf_var = 30;
ktype = 'rbf';

for m = 1:length(D)
    DD=[];
    Dm=D(m);
    for j = 1:length(D(m))
        DD = [DD Dm(j).D];
    end
    D(m).D=DD(:,randsample(size(DD,2),size(DD,1)));
end


d_Y = KMSC(D,Y,ktype,rbf_var);
cd(dirm);
save('d_Y.mat','d_Y');
cd(home);

cd(dirm);
open('d_Y.mat');
cd(home);
%% Documenting Results
gammas  = 0:1e-2:4.5;
gammas  = [gammas, 20];

% First 3 are alcyl2ft, alcycl3ft, alpipe
origT = t_Y;
t_Y(t_Y==1|t_Y==2|t_Y==3)= 0;
t_Y(t_Y~=0)= 1;

[gamk, Pcc, Pfa] = formROCS(d_Y,t_Y,gammas);


