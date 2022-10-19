input = linspace(-1,1,3);
i = 2;

N1 = 100;

lambda = [4 4; 4 4];
mu = [8 27; 5 21];
theta = [0.1 0.4; 0.1 0.3];
q = [15 17];
qD = [15 17; 15 17];
cost = [1 1];

[gOpt,r,~] = OptimalPolicySameEnvironments(120,120,lambda,mu,theta,10^(-input(i))*q,cost);
gOpt

%
WhitMat = IndexMatrixSame(mWhittleSame(120,lambda,mu,theta,10^(-input(i))*q,cost));
[gWhit] = PerfAnyPolSame(120,120,lambda,mu,theta,10^(-input(i))*q,cost,WhitMat);
%[gFixedEnvirD] = PerfAnyPolSame(120,120,lambda,mu,theta,10^(-input(i))*q,cost,IndexMatrixSame(mWhittle1EnvirLinearCost(120,mu,theta,cost)));
AverMat = IndexMatrixSame(mAverParamLinearCostSame(120,mu,theta,q,cost));
[gAverParam] = PerfAnyPolSame(120,120,lambda,mu,theta,10^(-input(i))*q,cost,AverMat);
    
hWhitAdd = (gWhit-gOpt) / gOpt * 100
%hFixedAddD = (gFixedEnvirD-gOpt) / gOpt * 100
hAverParam = (gAverParam-gOpt) / gOpt * 100
%}