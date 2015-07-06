%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hyperpsectral target detection using Kernel Matched Subspace
% Detectors; Final version for the manuscript submitted to IEEE 
% Trans. Image Processing.  
% Copyright (c) Heesung Kwon   11. 25, 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set what knid of energy term in the detection EQ. we use first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 5; % 1:k(y,y), 2: K(Ztb,kZtb)^2, 
            % 3:K(Ztt,Ztt)'ETET'K(Ztt,Ztt)+K(Zbb,Zbb)'EBEB'K(Zbb,Zbb)
            % 4:K(Ztb,Ztb)'*ETB*ETB'*K(Ztb,Ztb), ETB from a global background: KMSD
            % 5:Kernel matched filter: KTAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in a hyperspectral cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input = 3;
if input == 1
  B = 150; C = 129; R = 101;   %lower3
  fid = fopen('fr_test_sm_lower3.dat', 'r','b');
elseif input == 2
  B = 150; C = 574; R = 86;   %all targets
  fid = fopen('../Images/fr_all_targ.dat', 'r','b');
else
  %B = 180; C = 372; R = 135;   %drii 
  %fid = fopen('drii_ksubs.dat', 'r','b');
  B = 150; C = 372; R = 135;   %drii 
  fid = fopen('../Images/drii_bg_150.dat', 'r','b');
end
tmpdata = zeros(1, B*R*C);
tmpdata = fread(fid, 'ushort');
nc = max(tmpdata);
tmpdata = tmpdata./nc;
fclose(fid);
T = reshape(tmpdata,B,R,C);
clear tmpdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of eigen vectors for subspace estimation: AT, AB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input <= 2
  sst = 5; ssb = 10; %5  %forest: (5,5), drii:(2,6)best,(2,5),(1,6) 
else
  sst = 2; ssb = 6;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the target and background samples
% for subspace modeling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if input <= 2 % forest images (trees + grass)
  rbf_var = 10;  %10 used to be 30, 40 a possibility
  % Get first target samples
  NT1 = 18;
  fid = fopen('../Smpl_Spectra/samp_targ_18.dat','r');
  XT1 = fread(fid,'ushort');
  fclose(fid);
  XT1 = reshape(XT1, B, NT1);
  % Get second target samples
  %NT2 = 20;
  % read in target samples, the right most target
  %fid = fopen('samp_targ_20_rm.dat','r');
  %XT2 = fread(fid,'ushort');
  %fclose(fid);
  %XT2 = reshape(XT2, B, NT2);

  % Use only the first target
  NT = NT1;
  XT = [XT1];

  % Get background samples
  NB = 1656;  
  fid = fopen('../Smpl_Spectra/samp_bg_1656.dat','r');
  XB = fread(fid,'ushort');
  fclose(fid);
  XB = reshape(XB, B, NB);
else      % drii images (desert environments)
  rbf_var = 10;  % used to be 30, 40 a possibility
  % Get first target samples
  NT1 = 18;  %30; %30
  fid = fopen('../Smpl_Spectra/samp_targ_18_1fromleft_drii_150.dat','r');
  %fid = fopen('samp_targ_30_2fromleft_drii.dat','r');
  %fid = fopen('samp_targ_15_1fromleft_drii.dat','r');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply K-means to the background samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
km = 1;
if km == 1,
  %Select representative bg vectors using kmeans
  ncentres = NT*20; %3  % used to be 3,5 but 5 is sometimes better
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of the samples                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizet = size(XT,2); sizeb = size(XB,2);
sizetb = sizet+sizeb;
XT = XT/(nc); XB = XB/(nc);
XTB = [XT XB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate target template(matched filter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_sig = mean(XT,2);
XtB = [t_sig XB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Kernel matrix K(Xb,Xt), K(Xb,Xb), K(Xt,Xt), K(Xtb,Xtb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:sizeb,
for n=1:sizet,
  Kbt(m,n) = exp(-norm(XB(:,m)-XT(:,n))^2/rbf_var);
end
end
% K(Xb,Xb)
for m=1:sizeb,
for n=m:sizeb,
  Kbb(m,n) = exp(-norm(XB(:,m)-XB(:,n))^2/rbf_var);
  Kbb(n,m) = Kbb(m,n);
end
end
% K(Xt,Xt)
for m=1:sizet,
for n=1:sizet,
  Ktt(m,n) = exp(-norm(XT(:,m)-XT(:,n))^2/rbf_var);
end
end
for m=1:sizeb+sizet,
for n=m:sizeb+sizet,
  KTB(m,n) = exp(-norm(XTB(:,m)-XTB(:,n))^2/rbf_var);
  KTB(n,m) = KTB(m,n);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Kernel matrix KtBtB=K(XtB,XtB),KtBt=K(XtB,t_sig),KBt=K(XB,t_sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:sizeb+1,
for n=m:sizeb+1,
  KtBtB(m,n) = exp(-norm(XtB(:,m)-XtB(:,n))^2/rbf_var);
  KtBtB(n,m) = KtBtB(m,n);
end
end

for m=1:sizeb+1,
  KtBt(m) = exp(-norm(XtB(:,m)-t_sig)^2/rbf_var);
end
KtBt = KtBt';
for m=1:sizeb,
  KBt(m) = exp(-norm(XB(:,m)-t_sig)^2/rbf_var);
end
KBt = KBt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Centering of the kernel matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit_t = ones(sizet)/sizet;
unit_b = ones(sizeb)/sizeb;
unit_tb = ones(sizeb+sizet)/(sizeb+sizet);
unit_tB = ones(sizeb+1)/(sizeb+1);
Kbt = Kbt  - Kbt*unit_t - unit_b*Kbt + unit_b*Kbt*unit_t;
Ktb = Kbt';
Kbb = Kbb - Kbb*unit_b - unit_b*Kbb + unit_b*Kbb*unit_b;
Ktt = Ktt - Ktt*unit_t - unit_t*Ktt + unit_t*Ktt*unit_t;
KTB = KTB - KTB*unit_tb - unit_tb*KTB + unit_tb*KTB*unit_tb;
% matched filter
KtBtB = KtBtB - KtBtB*unit_tB - unit_tB*KtBtB + unit_tB*KtBtB*unit_tB; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA for the kernel matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EB, DB] = eig(Kbb);
[ET, DT] = eig(Ktt);
[ETB1, DTB] = eig(KTB);
% matched filter
[EtB, DtB] = eig(KtBtB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of the eigenvectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DT = real(diag(DT));
DB = real(diag(DB));
DTB = real(diag(DTB));
%matched filter
DtB = real(diag(DtB));

for m=1:sizet,
  ET(:,m) = ET(:,m)/(sqrt(DT(m)));
end
for m=1:sizeb,
  EB(:,m) = EB(:,m)/(sqrt(DB(m)));
end
for m=1:sizeb+sizet,
  ETB(:,m) = ETB1(:,m)./(sqrt(DTB(m)));
end
% matched filter
for m=1:sizeb+1,
  EtB(:,m) = EtB(:,m)./(sqrt(DtB(m)));
end

[at,bt] = sort(DT);
[ab,bb] = sort(DB);
[atb,btb] = sort(DTB);
% matched filter
[atB,btB] = sort(DtB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVs for separate energy normalization: ETT, EBB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bt(sizet) == sizet
  AT = ET(:,bt(sizet)-(sst-1):bt(sizet)); 
  ETT = ET(:,bt(sizet)-10:bt(sizet));  
else
  AT = ET(:,bt(sizet):bt(sizet)+sst-1);
  ETT = ET(:,bt(sizet):bt(sizet)+10);
end

if bb(sizeb) == sizeb
  AB = EB(:,bb(sizeb)-(ssb-1):bb(sizeb));
  EBB = EB(:,bb(sizeb)-10:bb(sizeb));
else
  AB = EB(:,bb(sizeb):bb(sizeb)+ssb-1);
  EBB = EB(:,bb(sizeb):bb(sizeb)+10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVs for joint energy normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jointrun = 15;
if btb(sizetb) == sizetb
  ETB_ = ETB(:,btb(sizetb)-jointrun:btb(sizetb));
else
  ETB_ = ETB(:,btb(sizetb):btb(sizetb)+jointrun);
end
% 8, 20 % 5, 20 %14,30 doesn't work
%3,3 also okay
% If use too small no of EV (w large EVal)
%  it causes problem
% 8, 10 better, %2,2 also okay
% drii
%ETB_ = ETB(:,1:10);

% matched filter
if btB(sizeb+1) == sizeb+1
  EtB_ = EtB(:,btB(sizeb+1)-jointrun:btB(sizeb+1));
else
  EtB_ = EtB(:,btB(sizeb+1):btB(sizeb+1)+jointrun);
end


%%%%%%%%%%%%%%%%%%%%%%%
% Mixed matrix creation 
%%%%%%%%%%%%%%%%%%%%%%%
mixed_mat = zeros(sst+ssb);
mixed_mat(1:sst,1:sst) = AT'*Ktt*AT;
mixed_mat(sst+1:(sst+ssb),sst+1:(sst+ssb)) = AB'*Kbb*AB;
mixed_mat(1:sst,sst+1:(sst+ssb)) = AT'*Ktb*AB;
mixed_mat(sst+1:(sst+ssb),1:sst) = AB'*Kbt*AT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel matched subspace detection main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enlarged hyperspectral cube for local background method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margin = 60;  % margin added to the whole cube 
T_local = zeros(B, R+margin, C+margin);
T_local(:, (margin/2):(margin/2+R-1), (margin/2):(margin/2+C-1)) = T;

if input == 1
  starty = 2; startx = 2; endy = 54; endx = 110;  %lower3
elseif input == 2
  starty = 1; startx = 1; endy = 78; endx = 460;  %endx = 564;  %all_targ
else
  starty = 20; startx = 50; endy = R-2; endx = C-2;  %drii
end

for i = starty:endy  
for j = startx:endx 

  y = T(:,i,j);
  if method <= 4
    % Calculate  KyXb
    for m=1:sizeb,
      KyXb(m,1) = exp(-norm(XB(:,m) - y)^2/rbf_var);
    end
    for m=1:sizet,
      KyXt(m,1) = exp(-norm(XT(:,m) - y)^2/rbf_var);
    end
    for m=1:sizet+sizeb,
      KyXtb(m,1) = exp(-norm(XTB(:,m) - y)^2/rbf_var);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KyXb = (KyXb - mean(KyXb));
    KyXt = (KyXt - mean(KyXt));
    KyXtb = (KyXtb - mean(KyXtb));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LT = [KyXt'*AT KyXb'*AB];
    RT = [AT'*KyXt; AB'*KyXb];
  else % matched filter
    for m=1:sizeb,
      KBy(m,1) = exp(-norm(XB(:,m) - y)^2/rbf_var); % KBy=KyXb;
    end
    for m=1:sizeb+1,
      KtBy(m,1) = exp(-norm(XtB(:,m) - y)^2/rbf_var); % KBy=KyXb;
    end
  end

  if method == 1  %k(y,y)
    energy = 1;
  elseif method == 2
    energy = KyXtb'*KyXtb;
  elseif method == 3
    energy = KyXt'*ETT*ETT'*KyXt + KyXb'*EBB*EBB'*KyXb;
  elseif method == 4
    energy = KyXtb'*ETB_*ETB_'*KyXtb;
  else % matched filter
    energynu = KtBt'*EtB_*EtB_'*KtBy;
    energyde = KtBt'*EtB_*EtB_'*KtBt;
  end

  if method <= 4
    term1 =  KyXb'*AB*AB'*KyXb;
    term2 = LT*pinv(mixed_mat)*RT;
    feature = (energy-term1)/(energy-term2);
    nu = (energy-term1);
    de =  (energy-term2);
  else % matched filter
    term1 = KBt'*AB*AB'*KBy;
    term2 = KBt'*AB*AB'*KBt;
    feature = (energynu-term1)/(energyde-term2);
    nu = (energynu-term1);
    de =  (energyde-term2);
  end

  if input <= 2
    abundance((i),(j)) = feature;
    ori((i),(j)) = y(40);
    top((i),(j)) = nu;
    bottom((i),(j)) = de;
    topproj((i),(j)) = term1;
    bottproj((i),(j)) = term2;
  else
    abundance((i-19),(j-49)) = feature;
    ori((i-19),(j-49)) = y(40);
    top((i-19),(j-49)) = nu;
    bottom((i-19),(j-49)) = de;
    topproj((i-19),(j-49)) = term1;
    bottproj((i-19),(j-49)) = term2;
  end
end
  i
end
figure
imagesc(abundance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROC curves generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 4
  [KPD,KPFA,Thits,False,Tpixs] = roc_gen_paper(abundance, input, 1);
  figure, plot(KPFA,KPD);
else
  [KPD_MF,KPFA_MF,Thits,False,Tpixs] = roc_gen_paper(abundance, input, 1);
  figure, plot(KPFA_MF,KPD_MF);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

