%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hyperpsectral target detection using Matched Subspace
% Detectors in the original input domain; Final version
% for the manuscript submitted to IEEE Trans. Image 
% Processing.
% Copyright (c) Heesung Kwon   11. 25, 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 2; % 1: subspace signal, 
            % 2: matched filter (target template)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in input hyperspectral cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input = 2;
if input == 1
  B = 150; C = 129; R = 101;   %lower3
  fid = fopen('../Images/fr_test_sm_lower3.dat', 'r','b');
elseif input == 2
  B = 150; C = 574; R = 86;   %all targets
  fid = fopen('../Images/fr_all_targ.dat', 'r','b');
else
  %B = 180; C = 372; R = 135;   %drii
  %fid = fopen('drii_ksubs.dat', 'r','b');
  B = 150; C = 372; R = 135;   %drii
  fid = fopen('../Images/drii_bg_150.dat', 'r','b');
end

%B = 150; C = 145; R = 41;   % 5tar
%B = 150; C = 254; R = 74;   %fr_pannels
%fid = fopen('fr_anomal_5tar.dat', 'r','b');
%fid = fopen('fr_pannels.dat', 'r','b');

tmpdata = zeros(1, B*R*C);
tmpdata = fread(fid, 'ushort');
nc = max(tmpdata);
tmpdata = tmpdata./nc;
fclose(fid);
T = reshape(tmpdata,B,R,C);
clear tmpdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in target and background samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input <= 2 % forest images
  NT = 18; 
  % read in target and background samples
  %fid = fopen('samp_targ_18_no_shadow.dat','r');
  fid = fopen('../Smpl_Spectra/samp_targ_18.dat','r');
  XT = fread(fid,'ushort');
  fclose(fid);
  XT = reshape(XT, B, NT);

  NB = 1656; %310,531,1196,420
  %fid = fopen('samp_bg_shadow_420.dat','r');
  %fid = fopen('samp_bg_shadow_1196.dat','r');
  fid = fopen('../Smpl_Spectra/samp_bg_1656.dat','r');
  %fid = fopen('samp_bg_310.dat','r');
  %fid = fopen('samp_bg_531.dat','r');
  XB = fread(fid,'ushort');
  fclose(fid);
  XB = reshape(XB, B, NB);
else      % drii images (desert environments)
  % Get first target samples
  NT1 = 18;  %30; %30
  fid = fopen('../Smpl_Spectra/samp_targ_18_1fromleft_drii_150.dat','r');
  %NT1 = 30;
  %fid = fopen('samp_targ_30_2fromleft_drii.dat','r');
  XT1 = fread(fid,'ushort');
  fclose(fid);
  XT1 = reshape(XT1, B, NT1);
                                                                                                     
  NT = NT1;
  XT = [XT1];
                                                                                                     
  % Get background samples
  %NB = 1932;
  %fid = fopen('samp_bg_1932_drii.dat','r');
  NB = 1300;
  fid = fopen('../Smpl_Spectra/samp_bg_1300_drii_150.dat','r');
  XB = fread(fid,'ushort');
  fclose(fid);
  XB = reshape(XB, B, NB);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply K-means to background samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
km = 1;
if km == 1,
  %Select representative bg vectors using kmeans
  ncentres = NT*20;   % used to be 2,5
  centres = zeros(ncentres, B);
  % Set up vector of options for kmeans trainer
  %options = foptions;
  options(1)  = 1;                % Prints out error values.
  options(5) = 1;
  options(14) = 60;               % Number of iterations.
  % Train the centres from the data
  [centres, options, post] = kmeans(centres, XB', options);
  XB = centres';
end

XT = XT/nc; XB = XB/nc;
XTB = [XT XB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target template (matched filter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_sig = mean(XT,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the target and background covariance
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cov_xt = cov(XT');
cov_xb = cov(XB');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input <= 2
  sst = 5; ssb = 10; %forest:(5,10) not very dependent on no of EVs,drii:(2,13) (NT*5)
else
  sst = 2; ssb = 13;
end
[EB, DB] = eig(cov_xb);
[ET, DT] = eig(cov_xt);
DT = diag(DT);
DB = diag(DB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sorting eigenvalues DT and DB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[at,bt] = sort(DT);
[ab,bb] = sort(DB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose significant eigenvectors for the target and background
% subspaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bt(B) == B
  AT = ET(:,bt(B)-(sst-1):bt(B));
else
  AT = ET(:,bt(B):bt(B)+sst-1);
end
if bb(B) == B
  AB = EB(:,bb(B)-(ssb-1):bb(B));
else
  AB = EB(:,bb(B):bb(B)+ssb-1);
end
%AB = EB(:,B-ssb+1:B);
%AT = ET(:,B-sst+1:B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenated target and background subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TB = [AT AB];

BB = AB*AB';
QQ = TB*pinv(TB'*TB)*TB';
I = eye(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input == 1
  starty = 2; startx = 2; endy = 54; endx = 110;  %lower3
elseif input == 2
  starty = 1; startx = 1; endy = 78; endx = 460; %endx = 564;  %all_targ
else
  starty = 20; startx = 50; endy = R-2; endx = C-2;  %drii
end
%starty = 8; startx = 8; endy = 35; endx = 138; %5tar
%starty = 8; startx = 8; endy = 67; endx = 246;  %fr_pannels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = starty:endy  % 10
for j = startx:endx % 100,110

  % input
  y = T(:,i,j);

  if method == 1               % signal subspace
    NU = y'*(I - BB)*y;
    DE = y'*(I - QQ)*y;
  else                         % matched filter
    NU = x_sig'*(I - BB)*y;
    DE = x_sig'*(I - BB)*x_sig;
  end

  features = NU/DE;

  if input <= 2
    conf((i),(j)) = features;
  else
    conf((i-19),(j-49)) = features;
  end

end
i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(conf);
%colormap(gray);
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROC curves generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 1
  [PD,PFA,Thits,False,Tpixs] = roc_gen_paper(conf, input, 2);
  figure, plot(PFA,PD);
else
  [PD_MF,PFA_MF,Thits,False,Tpixs] = roc_gen_paper(conf, input, 2);
  figure, plot(PFA_MF,PD_MF);
end



