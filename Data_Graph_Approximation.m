function [H,Qi,Ai,B,alpha,beta,obj] = Data_Graph_Approximation(X,num_nei,lambda1,lambda2,lambda3,Class)
 % X表示每一行表示一个数据，以行为单位
Num_view = length(X);
Num_fea = size(X{1},1);


Wi{Num_view} = [];
for i = 1:Num_view
    Wi{i} = constructW_PKN(X{i}', num_nei, 1);
end

Ai{Num_view} = [];
for i = 1:Num_view
    Ai{i} = ones(Num_fea,Class);
end

B = ones(Num_fea,Class);
H = ones(Num_fea,Class);
alpha = (1/Num_view)*ones(Num_view,1);
beta = (1/Num_view)*ones(Num_view,1);


Y1 = 0;
Y2 = zeros(Num_view,1);
Y3 = 0;
Y4 = zeros(Num_view,1);

f = alpha;
f0 = f;



mu = 0.1;
rho = 1.01;
max_mu = 10^5;

max_iter = 100;
for iter = 1:max_iter
%     iter
    % Qi
    Qi{Num_view} = [];
    for i = 1:Num_view
        [ui,~,vi] = svd((X{i}'*Ai{i}),'econ');
        Qi{i} = ui*vi';       
    end
    clear ui; clear vi;
    
     % Ai
     ri = zeros(Num_fea,Class);
     rii = ri;
     for i = 1:Num_view
         ri = ri + beta(i)*B;
         for j = 1:Num_view
             rii = rii+beta(i)*beta(j)*Ai{j};            
         end
         riii = rii-beta(i)*beta(j)*Ai{i};
         Ai{i} = (2*X{i}*Qi{i}+2*lambda2*ri-lambda2*riii)*pinv(2*Qi{i}'*Qi{i}+2*lambda2*beta(i)^2*eye(Class));
         Ai{i} = max(Ai{i},0);
     end
     clear ri;
     clear rii;
     clear riii;
     
     % B
     ei = zeros(Num_fea,Num_fea);
     eii = zeros(Num_fea,Class);
     for i = 1:Num_view
         ei = ei + alpha(i)*Wi{i};
         eii = eii + beta(i)*Ai{i};
     end
     B = (lambda1*ei*H+lambda2*eii+lambda3*H)*pinv(lambda1*H'*H+lambda2*eye(Class)+lambda3*eye(Class));
     B = max(B,0);
     clear eii;
     
     % H
     [ui,~,vi] = svd(2*lambda1*ei'*B+lambda3*B,'econ');
     H = ui*vi';
     clear ui; clear vi;
     
     % alpha
     
     PXi = [];
     for i = 1:Num_view
         PXi = [PXi Wi{i}(:)];
     end
     BBi = PXi'*PXi;
     ci = B*H';
     hbi = 2*PXi'*ci(:);
     alpha = pinv(2*BBi+mu*ones(Num_view,1)*ones(Num_view,1)'+mu*eye(Num_view))*(hbi+mu*ones(Num_view,1)*(1-Y1/mu)+mu*(f-Y2/mu)); 
     f = max(alpha+Y2/mu,0);
     clear BBi; clear hbi; clear PXi;
     
     
     % update beta 
     PXi = [];
     for i = 1:Num_view
         PXi = [PXi Ai{i}(:)];
     end 
     BBi = PXi'*PXi;
     hbi = 2*PXi'*B(:);
     beta = pinv(2*BBi+mu*eye(Num_view)+mu*ones(Num_view,1)*ones(Num_view,1)')*(hbi+mu*ones(Num_view,1)*(1-Y3/mu)+mu*(f0-Y4/mu));
     f0 = max(alpha+Y4/mu,0);
     clear BBi; clear hbi; clear PXi;
     
           
      Y1 = Y1+mu*(alpha'*ones(Num_view,1)-1);
      Y2 = Y2+mu*(alpha-f);
      Y3 = Y3+mu*(beta'*ones(Num_view,1)-1);
      Y4 = Y4+mu*(beta-f0);
      mu = min(rho*mu, max_mu);
      
      obj1 = 0;
      yi = zeros(Num_fea,Num_fea);
      yii = zeros(Num_fea,Class);
      for i = 1:Num_view
          obj1 = obj1+norm(X{i}-Ai{i}*Qi{i}','fro')^2;
          yi = yi+alpha(i)*Wi{i};
          yii = yii+beta(i)*Ai{i};
      end
      
      obj(iter) = obj1+lambda1*norm(yi-B*H','fro')^2+lambda2*norm(yii-B,'fro')^2+lambda3*norm(B-H,'fro')^2;
      if iter > 2
          if abs(obj(iter)-obj(iter-1)) < 10^-9;
              break
          end
      end  
end







