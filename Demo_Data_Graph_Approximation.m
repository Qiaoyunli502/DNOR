clear all;
clc,clear;
currentFolder = pwd;
warning off;
addpath(genpath('ClusteringMeasure'));
addpath(genpath('Datasets'));
addpath(genpath('Func'));
load MSRC;
% X = data;
% Y = (truelabel{1})';

Num_view = length(X);
XX{Num_view} = [];
for i= 1:Num_view
%     YY = X{i}';
    XX{i} = NormalizeFea(X{i},1);
%     XX{i} = NormalizeFea(X{i},1);
end
Class = length(unique(Y));

 lambda1 = [10^-8];
 lambda2 = [10^-8];
 lambda3 = [10^-10];
 num_nei = [15];
bestrate = 0;
for i = 1:length(lambda1)
    for j = 1:length(lambda2)
        for k = 1:length(lambda3)
            for ij = 1:length(num_nei)
            [H,Qi,Ai,B,alpha,beta,obj] = Data_Graph_Approximation(XX,num_nei(ij),lambda1(i),lambda2(j),lambda3(k),Class);
            B = NormalizeFea(B,1);
            P_label = kmeans(B,Class,'emptyaction','singleton','replicates',20,'display','off');
            result = Clustering8Measure(P_label, Y);
            rate = result(7)*100;
            
          
             if rate > bestrate
                 bestrate = rate;
                 bestlambda1 = lambda1(i);
                 bestlambda2 = lambda2(j);
                 bestlambda3 = lambda3(k); 
                 bestnum_nei = num_nei(ij);
             end
             
            end             
        end
    end
end

 
 


