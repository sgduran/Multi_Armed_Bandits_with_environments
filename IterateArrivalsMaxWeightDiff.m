input = linspace(1.5,3.75,10);

length_input=length(input);

lambda = [1 1; 1 1];
mu = [5 6; 8 7];
theta = [0 0; 0 0];
thetaBis = [0.0001 0.0001; 0.0001 0.0001];
qD = [3 2; 3 2];

%W = WIndices(mu,thetaBis,q);
%W

Optvec=[];
Whitvec=[];
MaxWeightvec=[];

hWhit=[];
hMaxWeight=[];


for i=1:length_input
 %For the performance of the optimal policies  
  
 [gOpt,r,VOpt] = OptimalPolicyDifferentEnvironmentsQuad(100,100,input(i)*lambda,mu,theta,qD,[1 1]);
 
 %{
 m = mWhittleDiffNoAbandonQuad(100,input(i)*lambda,mu,thetaBis,qD,[1 1]);
 MatrixWhittle = IndexDiffNoAbandon(m,mu);
 [gWhit] = PerfAnyPolDiffQuad(100,100,input(i)*lambda,mu,theta,qD,[1 1],MatrixWhittle);
 %}
 
 mMW = mMaxWeight(100,mu,[1 1]);
 MatrixMW = IndexMatrixDiff(mMW);
 [gMW] = PerfAnyPolDiffQuad(100,100,input(i)*lambda,mu,theta,qD,[1 1],MatrixMW);
    
 Optvec = [Optvec gOpt];
 %Whitvec = [Whitvec gWhit];
 MaxWeightvec = [MaxWeightvec gMW];
 
 %hWhitAdd = (gWhit-gOpt) / gOpt * 100;
 %hWhit = [hWhit hWhitAdd];

 hMaxWeightAdd = (gMW-gOpt) / gOpt * 100;
 hMaxWeight = [hMaxWeight hMaxWeightAdd];
  
end

%hWhit
hMaxWeight

hold off
h = figure;

%plot(input,hWhit)
%hold on
plot(input,hMaxWeight)

title('MaxWeight policy, scaling arrivals')
xlabel('\alpha, where arrivals are \lambda = \alpha*[1 1,1 1]')
ylabel('(gX - gOpt)/gOpt * 100')

xticks(input)
legend('Whittle Gap','Max Weight Gap')
%{
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'MaxWeightArrivalsDiff','-dpdf','-r0')
%}

%{
function [W] = WIndices(mu,theta,q)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+q(1)+q(2)) / (theta(1,i)*theta(2,i)+theta(2,i)*q(1)+theta(1,i)*q(2)));
    end
end

end
%}