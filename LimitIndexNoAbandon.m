% Parameters
alpha = 1;
lambda = alpha*[4 4; 4 4];
mu = alpha*[8 5; 27 21];
q = alpha*[15 17];
qD = alpha*[15 17; 15 17];
theta = [1 1; 1 1];
thetaBis = [0.001 0.001; 0.001 0.001];


fprintf('First we run the model with the regular indices \n')
fprintf('These indices should tend to infinity in environment 2 in both bandits, as theta goes to 0.\n');
for i = 1:3
	% Obtain indices
	m1 = mWhittleSame(130,lambda,mu,theta*(10^(-i)),q,[1 1]);
	% Print indices
	fprintf('Indices for first 30 states and theta = %f \n',theta[1,1]);
	m1Trun = m1(:,1:30,:)
end

% Obtain indices for our heuristics
fprintf('Now we obtain indices using the scaled limit in environment 2 \n'.)
fprintf('These indices should be %d and %d for envir 2, \n',mu[2,1],mu[2,2]);
fprintf('and a finite value for envir 1. \n');
m1NA = mWhittleSameNoAbandon(130,lambda,mu,thetaBis,q,[1 1]);
m1NATrun = m1(:,1:30,:)

