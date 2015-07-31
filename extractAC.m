function [ usableAsps, Bk ] = extractAC( filtered_run,eps,sigN,upperf,f_s,stopbins )
%% Input
% dir - directory to find Target Data
% filtered_run - matrix containing filtered run from Pond or TREX to be
% converted to a target strength AC
% eps - threshold for useful aspects for classification
% dim - array with the dimensions desired for the output (aspect res, freq
% res)
%% Output
% useableAsps - AC matrix for aspects with meaningful power (absolute sum
% greater than eps

% Bk -- AC matrix for Background aspects collected from this object run...
run = double(filtered_run);
len = size(run,2);

%ctr = findArcCtr(run,0.08*mean(mean(abs(run).^2)));

%% Begin Unpacking 
AC  = abs(fft(run,[],1));

% Windowing Useable spectral power
binsize = f_s/size(AC,1);
AC      = AC(1:ceil(upperf/binsize),:);

% Setting up smaller data structures for decimated form
aAC = zeros(size(AC,1),stopbins);
% Decimate along aspect if aperture is smaller than number of samples we
% have...
AC = AC(:,sum(abs(AC),1)>1);

if(stopbins<size(AC,2))
    for i = 1:size(AC,1)
    aAC(i,:) = resample(double(AC(i,:)),stopbins,size(AC,2));
    end
else
    aAC = AC;
end

wid = 60;

% Select enhanced aspects in the in center (the 'norm' ones all have this)
% for AC and lowest power observations for Bk

switch eps
   case 'ctr'
       ctr = floor(ctr*stopbins/size(run,2));
       AC = aAC(:,ctr-wid:ctr+wid-1);
   case 'pwr'
       [~, order] = sort(sum(abs(aAC).^2,1));
       order = fliplr(order);
       order = order(1:floor(length(order)*0.650));
       AC = aAC(:,order);
   case 'var'
       [~, order] = sort(var(normc(aAC),0,1));
       order = fliplr(order);
       order = order(1:floor(length(order)*0.1450));
       AC = aAC(:,order);
   otherwise
       % Selecting top aspects by hardlimiting power
            AC = aAC(:,sum(abs(aAC).^2,1)>eps);
end

mP = mean(sum(abs(aAC).^2,1));
Bk = aAC(:,sum(abs(aAC).^2,1)<1/20*mP);
% If we accidentally lost everything with power thresholding... 
if(size(AC,2)<floor(0.05*size(aAC,2)))
        [~, order] = sort(sum(abs(aAC).^2,1));
        order = fliplr(order);
        order = order(1:floor(length(order)*0.050));
        AC = aAC(:,order);
end

% For each aspect, decimate frequency bins to desired signal length
rAC = zeros(sigN,size(AC,2));
rBk = zeros(sigN,size(Bk,2));

for i = 1:size(AC,2)
    rAC(:,i) = abs(resample(AC(:,i),sigN,size(AC,1)))';
end
for i = 1:size(Bk,2)
    rBk(:,i) = abs(resample(Bk(:,i),sigN,size(Bk,1)))';
end
% Simulating Sparser stopping by taking every tenth aspect with a wobble of
% +-5 stops
rAC = reSort(rAC,12,4);
rBk = reSort(rBk,12,4);
usableAsps = normc(rAC);%20*log10(abs(normc(rAC)));%
Bk = normc(rBk);
end

function [ctrStop]=findArcCtr(run,thresh)
    ctrStop = size(run,2)/2;
    ctrCol = size(run,1)-10;
    for stop = 1:size(run,2)
        firstReturn = find(abs(run(:,stop)).^2 > thresh,1,'first');
        if(~isempty(firstReturn)&& firstReturn<=ctrCol)
               ctrStop = stop;
               ctrCol  = firstReturn;

        end
    end
end

% Simulate jumping with a wobble uniformly spaced in +- nbhdsize
function [AC]=reSort(rAC,jump,nbhdsz)
AC = zeros(size(rAC,1),size(rAC,2));
wobble = 0;
i = 0;
    while size(rAC,2)>= (mod(jump*(i)+wobble,size(rAC,2)-1)+1)
        AC(:,i+1) = rAC(:,(mod(jump*(i)+wobble,size(rAC,2)-1)+1));
        rAC(:,(mod(jump*(i)+wobble,size(rAC,2)-1)+1)) = [];
        wobble = round((rand(1)-0.5)*nbhdsz);
        i = i +1;
    end 
    AC(:,sum(abs(AC),1)<=1e-6)=[];
end




