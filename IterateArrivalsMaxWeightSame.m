input = linspace(0.25,3,12);

length_input=length(input);

lambda = [1 1; 1 1];
mu = [5 6; 8 7];
theta = [0 0; 0 0];
thetaBis = [0.0001 0.0001; 0.0001 0.0001];
q = [3 2];

W = WIndices(mu,thetaBis,q);
W

Optvec=[];
Whitvec=[];
MaxWeightvec=[];

hWhit=[];
hMaxWeight=[];


for i=1:length_input
 %For the performance of the optimal policies  
  
 [gOpt,r,VOpt] = OptimalPolicySameEnvironmentsQuad(80,80,input(i)*lambda,mu,theta,q,[1 1]);
 
 m = mWhittleSameNoAbandonQuad(80,input(i)*lambda,mu,thetaBis,q,[1 1]);
 MatrixWhittle = IndexSameNoAbandon(m,mu);
 [gWhit] = PerfAnyPolSameQuad(80,80,input(i)*lambda,mu,theta,q,[1 1],MatrixWhittle);
 
 mMW = mMaxWeight(80,mu,[1 1]);
 MatrixMW = IndexMatrixSame(mMW);
 [gMW] = PerfAnyPolSameQuad(80,80,input(i)*lambda,mu,theta,q,[1 1],MatrixMW);
    
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

title('MaxWeight policy, scaling arrivals')
xlabel('\alpha, where arrivals are \lambda = \alpha*[1 1,1 1]')
ylabel('(gX - gOpt)/gOpt * 100')

xticks(input)
legend('Whittle Gap','Max Weight Gap')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'MaxWeightArrivals','-dpdf','-r0')

function [W] = WIndices(mu,theta,q)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+q(1)+q(2)) / (theta(1,i)*theta(2,i)+theta(2,i)*q(1)+theta(1,i)*q(2)));
    end
end

end