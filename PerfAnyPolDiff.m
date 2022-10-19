% ----------------- Value iteration - Optimal ------------------------

function[gD,VD]=PerfAnyPolDiff(N1,N2,lambda,mu,theta,q,b,I)

disp('Value Iteration - Cost for any policy')

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

q112=q(1,1); % Rate for first environment going from state 1 to state 2
q121=q(1,2); % Rate for first environment going from state 2 to state 1
q212=q(2,1); % Rate for second environment going from state 1 to state 2
q221=q(2,2); % Rate for second environment going from state 2 to state 1


%For stability
if rho(lambda,mu,q) >= 1 
        
    %disp('No stable')
    %return
        
end

gamma = max(lambda11,lambda21) + max(lambda12,lambda22) + max(mu11,mu21) + max(mu12,mu22)...
    + N1 * max(theta11,theta21) + N2 * max(theta12,theta22) + max(q112,q121) + max(q212,q221); % Uniformization parameter

%Probabilities of going from one environment to the other in the time given
%by the exponential r.v. with rate gamma.


% coefficients for a linear cost function (b linear part)
b1 = b(1);
b2 = b(2);

precgamma = 10^(-8)*gamma; % numerical precision

VD = zeros(N1+1,N2+1,2,2); % Initializing value function
ConvergencePrecision = 10^-6; % The algorithm will stop after this precision is reached
MaxNumIterations = 80000;     %Maximum number of iterations if precision is not reached

iter = 0;

while (iter<MaxNumIterations)
    VDold = VD;
    VD = zeros(N1+1,N2+1,2,2);
    MinDiff = realmax;
    MaxDiff = -realmax;
        for j = 0:N1 % j is the state of bandit 1
            for k = 0:N2 % k is the state of bandit 2
                %revisar si esta bien que este el termino q(3-d)
                for d = 1:2
                for e = 1:2
                    VD(j+1,k+1,d,e) = gamma * cost(j,k,b);
                    
                        
                    VD(j+1,k+1,d,e) = VD(j+1,k+1,d,e) + q(1,d) * VDold(j+1,k+1,3-d,e)...
                        + q(2,e) * VDold(j+1,k+1,d,3-e)...
                        + lambda(d,1)*VDold(min(j+1,N1)+1,k+1,d,e)...
                        + lambda(e,2)*VDold(j+1,min(k+1,N2)+1,d,e)...
                        + j * theta(d,1)*VDold(max(j-1,0)+1,k+1,d,e)...
                        + k * theta(e,2)*VDold(j+1,max(k-1,0)+1,d,e);
                                                    
                    if I(j+1,k+1,d,e) == 2
                     VD(j+1,k+1,d,e) = VD(j+1,k+1,d,e) ...
                         + mu(e,2)*VDold(j+1,max(k-1,0)+1,d,e) ...
                         + (gamma-q(1,d)-q(2,e)-lambda(d,1)-lambda(e,2)- j*theta(d,1) ...
                         - k*theta(e,2)-mu(e,2)) * VDold(j+1,k+1,d,e);
                    elseif I(j+1,k+1,d,e) == 1
                     VD(j+1,k+1,d,e) = VD(j+1,k+1,d,e) ... 
                         + mu(d,1)*VDold(max(j-1,0)+1,k+1,d,e) ...
                         + (gamma-q(1,d)-q(2,e)-lambda(d,1)-lambda(e,2)-j*theta(d,1)...
                         - k*theta(e,2)-mu(d,1))*VDold(j+1,k+1,d,e);
                    end
                        
                    
                    end
                end
            end
        end
        
        
        VD = VD/gamma;
    % save VD preventively (for the case it gets stuck)
    if mod(iter,1000)==0
        %save VDiff VD;
    end

    if mod(iter,10)==0               %Every tenth iteration we update MinDiff and MaxDiff
        MinDiff = min(min(min(min(VD-VDold))));
        MaxDiff = max(max(max(max(VD-VDold))));
        if (MaxDiff-MinDiff)<ConvergencePrecision % algorithm has converged
            gD=(MaxDiff+MinDiff)/2;
            %save gDiff gD;
            
            fprintf('Value iteration converged in %3d iterations\n',iter+1);
            return
        end
    end
            
    if iter==(MaxNumIterations-1) % algorithm reached MaxNumIterations
        gD=(MaxDiff+MinDiff)/2;
        %save gDiff gD;
        %gOpt
        if MaxNumIterations >= 20000000
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
for a= 1:2
    M(a) = q(a,1) + q(a,2) ; 
    for b= 1:2
        phi(a,b) = q(a,3-b) / M(a) ;
    end


r = r + (phi(a,1) * lambda(1,a) + phi(a,2) * lambda (2,a)) / ...
    (phi(a,1) * mu(1,a) + phi(a,2) * mu(2,a)) ;

end

end
