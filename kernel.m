function kout = kernel(matx1, matx2, ktype, rbf_var)

kernel_types = {'rbf', 'imq', 'polynomial', 'sigmoidal', 'cosine'};

switch ktype

case 'rbf'
  %rbf_var = 30;
  kout = exp(-norm(matx1-matx2)^2/rbf_var);
case 'imq'
  kout = 1/sqrt(norm(matx1-matx2)^2+1);
case 'polynomial'
  theta = 1; const = 5;
  kout = (matx1'*matx2 + theta)^const;
case 'sigmoidal'
  theta1 = 1; theta2 = 1;
  kout = tanh(theta1*matx1'*matx2 + theta2);
case 'cosine'
  %kout = real(acos((matx1'*matx2)/(norm(matx1)*norm(matx2))));
  kout = (matx1'*matx2)/(norm(matx1)*norm(matx2));
otherwise
  error(['Unknown kernel type', ktype]);
end

% Centering kernel Matrix:
oneM = 1./size(kout,1)*ones(size(kout));
kout = (kout - oneM*kout - kout*oneM+ oneM*kout*oneM);
end

