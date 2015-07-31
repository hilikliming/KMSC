function [ gamk, P_cc, P_fa, m_Y ] = formROCS(d_Y,t_Y,gammas,dirm,tag,T,NT)
home = cd;
% Initializing min discriminant value vector for K observations and a
% decision vector m_P
K = length(t_Y);
jmin_T  = zeros(K,1);
jmin_NT = zeros(K,1);

m_Y = zeros(K,1);


% Finding minimal value from UXO and non UXO families of classes
for k = 1:K
% Strict Class decision
[~, m_Y(k)] = min(d_Y(k,:));
% Used in ratio comparison test for binary classifier
jmin_T(k) = min(d_Y(k,T));%/norm(d_YSVD(k,T));
jmin_NT(k) = min(d_Y(k,NT));%/norm(d_YSVD(k,NT));
end

% Setting UXO class I for objects classified as 1,2 or 3, other non UXO
% detections become non UXO class 0

%origM = m_Y;


for t = NT
t_Y(t_Y==t) = 0;
end


t_Y(t_Y~=0)   = 1;

%% ROC Curve Forming
P_cc     = zeros(length(gammas),1);
P_fa    = P_cc;
gam_i   = 1;
for gamma = gammas
    d  = zeros(K,1);
    for k = 1:K
        d(k)  = (jmin_T(k)/jmin_NT(k) < gamma);
    end
    d  = logical(d);
    P_cc(gam_i)   = sum(d & logical(t_Y))/sum(logical(t_Y));
    P_fa(gam_i)  = sum(d & ~logical(t_Y))/sum(~logical(t_Y));
    gam_i = gam_i +1;
end
[~,gamk]     = min(abs(1-(P_cc+P_fa)));
res  = [P_cc(gamk), P_fa(gamk)];


%% Displaying Result/Method Comparison and Decision Boundary
figure;
hold on
%       % NT stats         % Target Stats
plot(jmin_NT(t_Y==0),jmin_T(t_Y==0),'go'); %NT samples
%       % NT stats         % Target Stats
plot(jmin_NT(t_Y==1),jmin_T(t_Y==1),'rx'); %T Sample
plot((0:(gammas(2)-gammas(1)):1),gammas(gamk)*(0:(gammas(2)-gammas(1)):1)); %Decision Boundary
hold off
title([tag,' Decision Plane']);
legend('{NT}','{T}');%legend('J_1','J_2','J_3','J_4');%legend('J_1','J_2','J_3','J_4','J_C');%
axis([0,max(max(d_Y(:,NT))),0,max(max(d_Y(:,T)))]);


% Plotting ROCs and labeling Knee-Points
hndl=figure;
hold on
plot(P_fa,P_cc);%plot(P_faSVD,P_dSVD,P_faKSVD,P_dKSVD);%
legend(tag);%legend('SVD','K-SVD');%
tag1 = ['ROC for ',tag,' Tests'];
title(tag1); xlabel('P_{FA} (%)'); ylabel('P_{CC} (%)');
plot(P_fa(gamk), P_cc(gamk),'o');
axis([0, 1, 0, 1]);
hold off
cd(dirm);
%print(hndl,['TREXfil',num2str(sigA),'aMSDhd10m.png'],'-dpng');
savefig(hndl,['TREXfil',tag,'.fig']);
cd(home);


end

