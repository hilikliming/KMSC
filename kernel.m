function kout = kernel(matx1, matx2, ktype, rbf_var)

kernel_types = {'rbf', 'imq', 'polynomial', 'sigmoidal', 'cosine'};
kout = zeros(size(matx1,2),size(matx2,2));
switch ktype

case 'rbf'
  %rbf_var = 30;
  for i = 1:size(matx1,2)
      for j = 1:size(matx2,2)
            kout(i,j) = exp(-norm(matx1(:,i)-matx2(:,j))^2/rbf_var);
      end
  end
case 'imq'
    for i = 1:size(matx1,2)
      for j = 1:size(matx2,2)
            kout(i,j) = 1/sqrt(norm(matx1(:,i)-matx2(:,j))^2+1);
      end
    end
case 'polynomial'
    for i = 1:size(matx1,2)
      for j = 1:size(matx2,2)
        theta = 1; const = 5;
        kout(i,j) = (matx1(:,i)'*matx2(:,j) + theta)^const;
      end
    end
case 'sigmoidal'
    for i = 1:size(matx1,2)
      for j = 1:size(matx2,2)
        theta1 = 1; theta2 = 1;
        kout(i,j) = tanh(theta1*matx1(:,i)'*matx2(:,j) + theta2);
      end
    end
case 'cosine'
    for i = 1:size(matx1,2)
      for j = 1:size(matx2,2)
  %kout = real(acos((matx1'*matx2)/(norm(matx1)*norm(matx2))));
        kout(i,j) = (matx1(:,i)'*matx2(:,j))/(norm(matx1(:,i))*norm(matx2(:,j)));
      end
    end
otherwise
  error(['Unknown kernel type', ktype]);
end

% Centering kernel Matrix:
if size(kout,2)==size(kout,1)
    oneM = 1./size(kout,2).*ones(size(kout,2)); % kout in R^(MxN)
    kout = (kout - oneM*kout - kout*oneM+ oneM*kout*oneM);
else
    for vec = 1:size(kout,2)
        kout(:,vec) = kout(:,vec)-ones(size(kout(:,vec),1),1)*sum(kout(:,vec))/length(kout(:,vec));
    end
end

end

