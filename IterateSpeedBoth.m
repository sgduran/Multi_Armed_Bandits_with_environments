input = linspace(0.25,3,12);

length_input=length(input);

lambda = [6 6; 6 6];
mu = [14 32; 9 26];
theta = [0.1 0.1; 0.1 0.4];
q = [12 15];
qD = [12 15; 12 15];
% I LEAVE THIS WITH NORMAL COST
cost = [1 1];

W = WIndices(mu,theta,q);
W

Optvec=[];
Whitvec=[];
hWhit=[];

OptvecD=[];
WhitvecD=[];
hWhitD=[];


for i=1:length_input
 %For the performance of the optimal policies  
  
[gOpt,r,VOpt] = OptimalPolicySameEnvironments(100,100,lambda,mu,theta,input(i)*q,cost);
[gWhit] = PerfAnyPolSame(100,100,lambda,mu,theta,input(i)*q,cost,IndexMatrixSame(mWhittleSame(100,lambda,mu,theta,input(i)*q,cost)));
      
Optvec = [Optvec gOpt];
Whitvec = [Whitvec gWhit];
 
hWhitAdd = (gWhit-gOpt) / gOpt * 100;

hWhit = [hWhit hWhitAdd];
 
[gOptD,rD,VOptD] = OptimalPolicyDifferentEnvironments(100,100,lambda,mu,theta,input(i)*qD,cost);
 
[gWhitD] = PerfAnyPolDiff(100,100,lambda,mu,theta,input(i)*qD,cost,IndexMatrixDiff(mWhittleDiff(100,lambda,mu,theta,input(i)*qD,cost)));
    
OptvecD = [OptvecD gOptD];
WhitvecD = [WhitvecD gWhitD];
 
hWhitAddD = (gWhitD-gOptD) / gOptD * 100;
hWhitD = [hWhitD hWhitAddD];
 
end

hWhit
hWhitD

hold off
h = figure;

plot(input,hWhit)
hold on
plot(input,hWhitD)

title('Optimal gap, scaling speed')
xlabel('\alpha, where speed is q = alpha*[12 15,12 15]')
ylabel('(gX - gOpt)/gOpt * 100')

xticks(input)
legend('Whittle Gap','Whittle Diff Gap')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'SpeedBothEx9','-dpdf','-r0')

function [W] = WIndices(mu,theta,q)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+q(1)+q(2)) / (theta(1,i)*theta(2,i)+theta(2,i)*q(1)+theta(1,i)*q(2)));
    end
end

end
