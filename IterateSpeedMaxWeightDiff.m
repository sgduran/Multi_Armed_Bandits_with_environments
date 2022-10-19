input = linspace(0.25,3,12);

length_input=length(input);

lambda = 1.75*[1 1; 1 1];
mu = [5 6; 8 7];
theta = [0 0; 0 0];
thetaBis = [0.0001 0.0001; 0.0001 0.0001];
qD = [3 2; 3 2];

W = WIndices(mu,theta,qD);
W

Optvec=[];
Whitvec=[];
MaxWeightvec=[];

hWhit=[];
hMaxWeight=[];


for i=1:length_input
 %For the performance of the optimal policies  
  
 [gOpt,r,VOpt] = OptimalPolicyDifferentEnvironmentsQuad(80,80,lambda,mu,theta,input(i)*qD,[1 1]);
 
 m = mWhittleDiffNoAbandonQuad(80,lambda,mu,thetaBis,input(i)*qD,[1 1]);
 MatrixWhittle = IndexDiffNoAbandon(m,mu);
 [gWhit] = PerfAnyPolDiffQuad(80,80,lambda,mu,theta,input(i)*qD,[1 1],MatrixWhittle);
 
 mMW = mMaxWeight(80,mu,[1 1]);
 MatrixMW = IndexMatrixDiff(mMW);
 [gMW] = PerfAnyPolDiffQuad(80,80,lambda,mu,theta,input(i)*qD,[1 1],MatrixMW);
    
 Optvec = [Optvec gOpt];
 Whitvec = [Whitvec gWhit];
 MaxWeightvec = [MaxWeightvec gMW];
 
 hWhitAdd = (gWhit-gOpt) / gOpt * 100;
 hWhit = [hWhit hWhitAdd];

 hMaxWeightAdd = (gMW-gOpt) / gOpt * 100;
 hMaxWeight = [hMaxWeight hMaxWeightAdd];
  
end

hWhit
hMaxWeight

hold off
h = figure;

plot(input,hWhit)
hold on
plot(input,hMaxWeight)

title('MaxWeight policy, scaling speed. Separate Envirs')
xlabel('\alpha, where speed is q = \alpha*[3 2]')
ylabel('(gX - gOpt)/gOpt * 100')

xticks(input)
legend('Whittle Gap','Max Weight Gap')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'MaxWeightSpeedDiff','-dpdf','-r0')

function [W] = WIndices(mu,theta,q)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+q(1)+q(2)) / (theta(1,i)*theta(2,i)+theta(2,i)*q(1)+theta(1,i)*q(2)));
    end
end

end