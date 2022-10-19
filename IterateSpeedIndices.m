%WE TAKE SPEED TO 0 (AND SOMETIMES THETA TO 0) TO CHECK IF WE CAN DERIVE
%EXPRESSIONS FOR THE INDICES.

input = [1 0.01 0.0001];
%input = [0.000001];
%input = [0.001];


%{
lambda = [3 3; 3 3];
mu = [1 10; 4.8 9.5];
theta = [0.1 0.1 ; 0.3 0.3];
q = [2 4];
%}

lambda = [2 2; 2 2];
mu = [6 10; 10 6];
%theta = [0.00001 0.00001; 0.00001 0.00001];
theta = [1 1; 1 1];
theta = 0.001 * theta;
q = [2 2];

length_input=length(input);


for i=1:length_input

 
 m = mWhittleSame(80,lambda,mu,theta,input(i)*q,[1 1]);
%CLASSIC INDICES
 m(:,1:15,1)
 m(:,1:15,2)
 
 m = mWhittleSameNoAbandon(80,lambda,mu,theta,input(i)*q,[1 1]);
 %INDICES MULTIPLIED BY THETA WHEN GOOD ENVIR
 mNA = m(:,1:15,1)
 mNA = m(:,1:15,2)
end
disp end