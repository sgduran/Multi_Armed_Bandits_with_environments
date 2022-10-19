%%% THIS IS TO PLOT THE INDICES WHEN THE SPEED GOES TO 0

%
input = [1 0.1 0.01 0.001];

lambda = [4 4; 4 4];
mu = [8 27; 5 21];
theta = [0.1 0.4; 0.1 0.3];
q = [15 17];

length_input=length(input);

mTrun=zeros(length_input,2,30,2);


for i=1:length_input
 
 m = mWhittleSame(80,lambda,mu,theta,input(i)*q,[1 1]);
 mTrun(i,:,:,:) = m(:,1:30,:);
 %Usar comando squeeze(mTrun(i,:,:,:)) para ver bien que queda
 
  
  
end
disp end
%}

hold off
h = figure;

index1 = zeros(2,2,length_input);
xaxis = -log10(input);
for k=1:2
    for d=1:2
        %plot(xaxis(1):xaxis(length(xaxis)),ones(xaxis(length(xaxis))-xaxis(1)+1) * mu(d,k)/theta(d,k))
        %hold on
        %EVOLUTION OF THE INDEX OF STATE N = 1
        state = 1;
        index1(d,k,:) = squeeze(mTrun(:,k,state+1,d));
        %plot(xaxis,index1);
        %hold on
    end
end

h = figure;

plot(xaxis,squeeze(index1(1,1,1:length_input)),xaxis,squeeze(index1(2,1,1:length_input)),':',xaxis,squeeze(index1(1,2,1:length_input)),'--',xaxis,squeeze(index1(2,2,1:length_input)),'-.');

title(sprintf('Covergence of the index for state m=%d',state),'FontSize',14)
xlabel('\beta', 'FontSize',22)
ylabel(sprintf('W_k(%d,d)',state),'FontSize',14)

%VECTOR OF LIMITS
limit = [mu(1,1)/theta(1,1) mu(2,1)/theta(2,1) mu(1,2)/theta(1,2) mu(2,2)/theta(2,2)];
limitsort = sort(limit);
yticks(limitsort)
%set ( gca, 'xdir', 'reverse' )
xticks([0:length_input])
%yticklabels({sprintf('mu_k^{(d)} = %.2f',limitsort(1)),sprintf('mu_k^{(d)} = %.2f',limitsort(2)),sprintf('mu_k^{(d)} = %.2f',limitsort(3)),sprintf('mu_k^{(d)} = %.2f',limitsort(4))})
legend('k = 1, d = 1','k = 1, d = 2','k = 2, d = 1','k = 2, d = 2','Location','southeast','FontSize',10)
%plot(1:10);
%
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+0.1, pos(4)+0.5])
print(h,'SlowSpeed','-dpdf','-r0')
%}