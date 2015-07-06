function [pd,pfa,tothit,false,targpixs] = roc_gen(abun, input, algo);

NEG = -100;
if input <= 2  % forest images
  % 14 target locations 
  LY = [32; 30; 29; 30; 28; 30; 28; 32; 30; 29; 27; 28; 27; 29]*2;        % first ori (31,12) (33,16)
  LX = [13; 28; 43; 57; 70; 84; 99; 113; 141; 155; 168; 184; 200; 212]*2; % second ori(30,27) (33,32)
  RY = [33; 32; 33; 33; 32; 31; 30; 32; 31; 30; 30; 30; 30; 30]*2;
  RX = [16; 31; 45; 59; 71; 87; 100; 114; 142; 157; 169; 187; 201; 214]*2;
  num_targ = 14;
else
% LY = [34; 34; 34; 35; 36; 38]*2;
% LX = [41; 67; 95; 118; 141; 168]*2;
% RY = [35; 36; 38; 36; 39; 40]*2;
% RX = [45; 71; 97; 120;143; 171]*2;
  LY = [48; 49; 47; 50; 52; 59];
  LX = [33; 85; 141; 184; 233; 287];
  RY = [51; 53; 55; 53; 59; 62];
  RX = [39; 89; 143; 189; 236; 291];
  num_targ = 6;
% LY = LY-20; LX=LX-50; RY=RY-20; RX=RX-50;
end

% No of total target pixels
targpixs = 0;
tmp = zeros(1, num_targ);
for i = 1:num_targ
  if input == 3;
    tmp(i) = (RX(i)-LX(i)+1)*(RY(i)-LY(i));
  else
    tmp(i) = (RX(i)-LX(i)+1)*(RY(i)-LY(i)+1);
  end
  targpixs = tmp(i) + targpixs;
  fprintf('%d target pixels: %d\n', i, tmp(i));
end

% target bin
targ_bin = zeros(num_targ,1);

% No of total pixels in abundance
[r, c] = size(abun);
totpixs = r*c;

% Normalization
abun = abun/max(max(abun));

for i = 1:3500 %1500
  [maxval(i), indx(i)] = max(max(abun));
  [maxval(i), indy(i)] = max(max(abun'));
  abun(indy(i), indx(i)) = NEG;
end

if algo == 1   % kernel
  th = maxval(30);  %60
  step = .005;   %.005
else
  th = maxval(30); %10
  step = .01;   %.2
end

%for i = 1:400  %200
%tothit = 0; false = 0;
%for j = 1:450 %150 

for i = 1:3000  %200
tothit = zeros(1,num_targ); false = 0;
for j = 1:3000 %150

  hit = 0;
  if maxval(j) >= th-(i-1)*step
    for m = 1:num_targ
       if (indx(j) >= (LX(m)-1) & indx(j) <= (RX(m)+1) & ... % -1 , +1
              indy(j) >= (LY(m)-1) & indy(j) <= (RY(m)+1))   % -1 , 0 also okay 
         hit = 1;
         if i == 200
           targ_bin(m) = targ_bin(m)+1;
         end
         tothit(m) = tothit(m) + 1;
         break;
       end 
    end
    if hit == 0
      false = false + 1;
    end
  else 
    break;
  end
end
  totnumhit = 0;
  for m=1:num_targ
    if tothit(m) > tmp(m)
       tothit(m) = tmp(m);
    end
  end
  for m=1:num_targ
    totnumhit = totnumhit + tothit(m);
  end
  pd(i) = totnumhit/targpixs;
  pfa(i) = false/totpixs;
end

for i=1:num_targ
  fprintf('%d target bin = %d\n', i, targ_bin(i));
end

