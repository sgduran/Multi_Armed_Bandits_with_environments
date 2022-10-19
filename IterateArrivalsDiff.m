%Small scale:
input = linspace(0.25,3,12);
%input = [input 1.5 2];
%Large scale:
%input=[1 10 100 1000];

length_input=length(input);

lambda = [4 4; 4 4];
mu = [8 27; 5 21];
theta = [0.1 0.4; 0.1 0.3];
qD = [15 17; 10 2];
cost = [1 1];

WD = WIndicesD(mu,theta,qD);
WD

OptvecD=[];
WhitvecD=[];
hWhitD=[];


for i=1:length_input
 %For the performance of the optimal policies  
  
 [gOptD,rD,VOptD] = OptimalPolicyDifferentEnvironments(120,120,input(i)*lambda,mu,theta,qD,cost);
 
 [gWhitD] = PerfAnyPolDiff(120,120,input(i)*lambda,mu,theta,qD,cost,IndexMatrixDiff(mWhittleDiff(120,input(i)*lambda,mu,theta,qD,cost)));
    
 OptvecD = [OptvecD gOptD];
 WhitvecD = [WhitvecD gWhitD];
 
%hWhitAdd = (Whitvec(i)-Optvec(i)) / Optvec(i) * 100;
hWhitAdd = (gWhitD-gOptD) / gOptD * 100;

hWhitD = [hWhitD hWhitAdd];
 
 %save Optvec1 Optvec1
 %save Whitvec Whitvec
 
 
  
end

OptvecD
WhitvecD
hWhitD

hold off
h = figure;

plot(input,hWhitD)

title('Optimality gap, scaling arrivals')
xlabel('\gamma','FontSize',18)
ylabel('(gX - gOpt)/gOpt * 100')

xticks(input)
legend('Whittle Opt Gap')

%
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'ArrivalsIndSameRates','-dpdf','-r0')
%}

function [W] = WIndicesD(mu,theta,qD)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+qD(i,1)+qD(i,2)) / (theta(1,i)*theta(2,i)+theta(2,i)*qD(i,1)+theta(1,i)*qD(i,2)));
    end
end

end