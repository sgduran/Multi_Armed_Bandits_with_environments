repeat = 1;
%for AB =1:repeat
W3 = 1
W2 = 0
while (W3 >= W2-0.001 || W3 < 0) && repeat < 1000

%
%CASE 1: OUR CONDITION;
r = randi([1 20],1,9);
opt = randi([1 1],1,2);
%
%r(3) = randi([1 10^opt(1)],1,1);
critmu2 = r(3) * (r(6)+r(7)+r(8)) / (r(5)+r(7)+r(8));
%r(4) = randi([ceil(critmu2) ceil(critmu2)+10^opt(2)],1,1);

r(4) = randi([ceil(max(critmu2,r(3))) max(ceil(min(critmu2,r(3)))+1,ceil(max(critmu2,r(3))))],1,1);

N1 = 100;
%
lambda = [r(1) r(2)];
mu = [r(3) r(4)];
theta = [r(5) r(6)];
q = [r(7) r(8)];

gamma = max(lambda(1),lambda(2)) + max(mu(1),mu(2)) + N1 * max(theta(1),theta(2)) + q(1) + q(2);


p=1;
W1 = (mu(p)) * ((theta(3-p)+q(1)+q(2)) / (theta(1)*theta(2)+theta(2)*q(1)+theta(1)*q(2)));
p=2;
W2 = (mu(p)) * ((theta(3-p)+q(1)+q(2)) / (theta(1)*theta(2)+theta(2)*q(1)+theta(1)*q(2)));

%rand : numero aleatorio entre 0 y 1
W = rand*W2;
%W = 0.5*W2;

iter = 100000;
C = [(W*(theta(2)+q(2))-mu(2)) / (q(2)*mu(2)) W/mu(2)];
[Vold,G,iter,t1] = PerfIndexStopInConv(N1,lambda,mu,theta,q,1,W,C,iter);

mu(2) / (theta(2)+gamma);
C
%(W*(theta(2)+q(2))-mu(2)) / (q(2)*mu(2))
G(1:2,1,iter)
iter
t1

W3 = (mu(1)*mu(2)) / (mu(1)*theta(2)+q(2)*(mu(1)-mu(2)))

steps = iter*0.1;
%USE THE FOLLOWING LINE IN TERMINAL (NOT HERE) TO ADD STEPS AFTER
%CONVERGENCE
%[Vold,G,iter] = PerfIndexSteps(N1,lambda,mu,theta,q,1,W,C,100000,Vold,G,iter,steps);
%G(1:2,1,iter)


GTiempo = squeeze(G(1,2,:));

if AB == repeat
    disp('Ok')
end

[C,I] = min(min(mu(2)*G(1:10,2,1:iter)-mu(1)*G(1:10,1,1:iter)));

repeat = repeat +1

end