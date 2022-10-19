% ----------------- Value iteration - Optimal ------------------------

function[g,V]=PerfAnyPolSame(N1,N2,lambda,mu,theta,q,b,I)

disp('Value Iteration - Cost for any policy')

%I : matrix which indicates where to serve in each (j,k)

mu11 = mu(1,1); % When state of env. is 1, departure rate for bandit 1
mu12 = mu(1,2); % When state of env. is 1, departure rate for bandit 2
mu21 = mu(2,1); % When state of env. is 2, departure rate for bandit 1
mu22 = mu(2,2); % When state of env. is 2, departure rate for bandit 2

lambda11 = lambda(1,1); % When state of env. is 1, arrival rate for bandit 1
lambda12 = lambda(1,2); % When state of env. is 1, arrival rate for bandit 2
lambda21 = lambda(2,1); % When state of env. is 2, arrival rate for bandit 1
lambda22 = lambda(2,2); % When state of env. is 2, arrival rate for bandit 2

theta11 = theta(1,1); % When state of env. is 1, abandon. rate for bandit 1
theta12 = theta(1,2); % When state of env. is 1, abandon. rate for bandit 2
theta21 = theta(2,1); % When state of env. is 2, abandon. rate for bandit 1
theta22 = theta(2,2); % When state of env. is 2, abandon. rate for bandit 2

q12=q(1); % Rate for environment going from state 1 to state 2
q21=q(2); % Rate for environment going from state 2 to state 1

%{
%For stability
if rho(lambda,mu,q) >= 1 
        
    disp('No stable')
    return
        
end
%}

gamma = max(lambda11,lambda21) + max(lambda12,lambda22) + max(mu11,mu21) + max(mu12,mu22)...
    + N1 * max(theta11,theta21) + N2 * max(theta12,theta22) + q12 + q21 ; % Uniformization parameter

%Probabilities of going from one environment to the other in the time given
%by the exponential r.v. with rate gamma.


% coefficients for a linear cost function (b linear part)
b1 = b(1);
b2 = b(2);

precgamma = 10^(-8)*gamma; % numerical precision

V = zeros(N1+1,N2+1,2); % Initializing value function
ConvergencePrecision = 10^-6; % The algorithm will stop after this precision is reached
MaxNumIterations = 80000;     %Maximum number of iterations if precision is not reached

iter = 0;

while (iter<MaxNumIterations)
    Vold = V;
    V = zeros(N1+1,N2+1,2);
    MinDiff = realmax;
    MaxDiff = -realmax;
        for j = 0:N1 % j is the state of bandit 1
            for k = 0:N2 % k is the state of bandit 2
                %revisar si esta bien que este el termino q(3-d)
                for d = 1:2
                    V(j+1,k+1,d) = gamma * cost(j,k,b);
                    
                    V(j+1,k+1,d) = V(j+1,k+1,d) + q(d) * Vold(j+1,k+1,3-d) ...
                        + lambda(d,1)*Vold(min(j+1,N1)+1,k+1,d)...
                        + lambda(d,2)*Vold(j+1,min(k+1,N2)+1,d)...
                        + j * theta(d,1)*Vold(max(j-1,0)+1,k+1,d)...
                        + k * theta(d,2)*Vold(j+1,max(k-1,0)+1,d);
                                                    
                    if I(j+1,k+1,d) == 2
                     V(j+1,k+1,d) = V(j+1,k+1,d) + mu(d,2)*Vold(j+1,max(k-1,0)+1,d) ...
                         + (gamma-q(d)-lambda(d,1)-lambda(d,2)-j*theta(d,1)...
                         - k*theta(d,2)-mu(d,2))*Vold(j+1,k+1,d);
                    elseif I(j+1,k+1,d) == 1
                     V(j+1,k+1,d) = V(j+1,k+1,d) + mu(d,1)*Vold(max(j-1,0)+1,k+1,d) ...
                         + (gamma-q(d)-lambda(d,1)-lambda(d,2)-j*theta(d,1)...
                         - k*theta(d,2)-mu(d,1))*Vold(j+1,k+1,d);
                    end
                    
                    
                end
                end
        end
        
        V = V/gamma;
    % save V preventively (for the case it gets stuck)
    if mod(iter,1000)==0
        %save VWhittleIndex V;
    end

    if mod(iter,10)==0               %Every tenth iteration we update MinDiff and MaxDiff
        MinDiff = min(min(min(V-Vold)));
        MaxDiff = max(max(max(V-Vold)));
        if (MaxDiff-MinDiff)<ConvergencePrecision % algorithm has converged
            g=(MaxDiff+MinDiff)/2;
            save g g;
            
            fprintf('Value iteration converged in %3d iterations\n',iter+1);
            return
        end
    end
            
    if iter==(MaxNumIterations-1) % algorithm reached MaxNumIterations
        g=(MaxDiff+MinDiff)/2;
        %save g g;
        %gOpt
        if MaxNumIterations >= 200000000
            display('Does not converge!')
            return
        end
            
        antw = 20000;
        if (antw >0 )
            MaxNumIterations=MaxNumIterations+antw;
        else 
            return
        end
    end
    
    iter=iter+1;
        
end
end
function c = cost(j,k,b)
% We have linear cost 
    c = (j)*b(1) + (k)*b(2);
    
end
function r = rho(lambda,mu,q)
%For the stability.

%For the stationary measure of the environments
%phi(a,b) is the stationary measure of the environment of class a in state b
r = 0;
phi = zeros(2,1);

M = q(1) + q(2) ; 
    for z= 1:2
        phi(z) = q(3-z) / M ;
    end

for a= 1:2
r = r + (phi(1) * lambda(1,a) + phi(2) * lambda (2,a)) / ...
    (phi(1) * mu(1,a) + phi(2) * mu(2,a)) ;
end

end
