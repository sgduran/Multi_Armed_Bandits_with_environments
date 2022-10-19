function[m]=mWhittleSame(N,lambdalong,mulong,thetalong,q,blong)

m = zeros(2,N+1,2);



for i=1:2
    lambda = lambdalong(:,i);
    mu = mulong(:,i);
    theta = thetalong(:,i);
    b = blong(i);
    [~,W,~,~,~,~] = WhittleIndexForcedBorders(N+50,lambda,mu,q,theta,b);
    for l=0:N
        for d=1:2
            m(i,l+1,d)= W(l+1,d);
            %WhittleIndex(l,N+50,lambda,mu,b,i,d);
        end
    end
end

end


function [What,W,I,Pol,A,k] = WhittleIndexForcedBorders(X0,lambda,mu,q,theta,b)
%X0 is where we truncate, b is value for the linear cost

Wbar = zeros(X0+2,X0+2,X0+2,X0+2);
A = zeros(3,3);
%What = zeros(1,2*X0);
What = [-1];
I = [2];
I_row = [-1];
I_col = [-1];
Index = [];
W = zeros(X0,2);

%for k=1:2*X0-1
k=1;
while I(k) > 1
    I(1)=0;
    for i=(I_row(k)):min((I_row(k)+2),X0-1)
        for j=(I_col(k)):min((I_col(k)+2),X0-1)
            Wbar(I_row(k)+2,I_col(k)+2,i+2,j+2) = LinearFunctions([I_row(k) I_col(k)],[i j],X0+50,lambda,mu,q,theta,b);
            if Wbar(I_row(k)+2,I_col(k)+2,i+2,j+2) < What(k)-(10^(-2))
            A(i-I_row(k)+1,j-I_col(k)+1) = Inf;
            else
            A(i-I_row(k)+1,j-I_col(k)+1) = Wbar(I_row(k)+2,I_col(k)+2,i+2,j+2);
            end
            %maybe A(3,x) or A(x,3) is not defined in the last for
        end
    end
    A(1,1) = Inf;
        
    if I_row(k) == X0-2
    %in case row is in the border, we have to discard possible minimums
        for j=1:3
            A(3,j) = Inf;
        end
    elseif I_row(k) == X0-1
        for i=2:3
            for j=1:3
                A(i,j) = Inf;
            end
        end
    end
    if I_col(k) == X0-2
    %in case column is in the border, we have to discard possible minimums
        for i=1:3
            A(i,3) = Inf;
        end
    elseif I_col(k) == X0-1
        for j=2:3
            for i=1:3
                A(i,j) = Inf;
            end
        end
    end
    
    [What(k+1),I(k+1)] = min(A(:));
    [I_row(k+1),I_col(k+1)] = ind2sub(size(A),I(k+1));
    I_row(k+1) = I_row(k+1) + I_row(k) - 1;
    I_col(k+1) = I_col(k+1) + I_col(k) - 1;
    
    %{
    for x = I_row(k):I_row(k+1)-1
        for y = I_col(k):X0
            if slope([x y],X0+50,lambda,mu,q,theta) < slope([I_row(k+1) I_col(k+1)],X0+50,lambda,mu,q,theta)
                %disp(['Condition is not working in iteration ',num2str(k) ...
                 %   ', ' num2str([I_row(k+1),I_col(k+1)]) ' is less steep than ', num2str([x,y]), '.'])
            end
        end
    end
    for y = I_col(k):I_col(k+1)-1
        for x = I_row(k):X0
            if slope([x y],X0+50,lambda,mu,q,theta) < slope([I_row(k+1) I_col(k+1)],X0+50,lambda,mu,q,theta)
                %disp(['Condition is not working in iteration ',num2str(k) ...
                %    ', ' num2str([I_row(k+1),I_col(k+1)]) ' is less steep than ', num2str([x,y]), '.'])
            end
        end
    end
    %}
    k = k + 1;

    
end
k = k-1;
What(k+1) = [];

for m=1:k-1
    for l=I_row(m)+1:I_row(m+1)
        W(l+1,1) = What(m+1);
    end
    for p=I_col(m)+1:I_col(m+1)
        W(p+1,2) = What(m+1);
    end
end

for i=1:k
    Pol(1,i)=I_row(i);
    Pol(2,i)=I_col(i);
end

%{
OLD: I force borders to be equal to the indices
Crit = [(mu(1)) * ((theta(2)+q(1)+q(2)) / (theta(1)*theta(2)+theta(2)*q(1)+theta(1)*q(2)))  (mu(2)) * ((theta(1)+q(1)+q(2)) / (theta(1)*theta(2)+theta(2)*q(1)+theta(1)*q(2)))];

for l=1:X0-1
    if W(l+1,1) > Crit(1) || W(l+1,1) == 0
        W(l+1,1) = Crit(1);
    end
end
for p=1:X0-1
    if W(p+1,2) > Crit(2) || W(p+1,2) == 0
        W(p+1,2) = Crit(2);
    end 
end
%}

%I force it just not to be 0
for l=5:X0-1
    if W(l+1,1) == 0
        W(l+1,1) = W(l,1);
    end
end
for p=5:X0-1
    if W(p+1,2) == 0
        W(p+1,2) = W(p,2);
    end 
end

end
function L = LinearFunctions(n,l,X1,lambda,mu,q,theta,b)
%X1 is where we truncate, 
%b is the vector for the linear cost

pi = StationaryMeasure(n,X1,lambda,mu,q,theta);

num = 0;
for j=0:X1
    for d=1:2
        num = num + pi(2*j+d) * cost2(j,b) ;
    end
end

den = 0;
m = min(n);
[M,I] = max(n);

for j=0:(2*m+1)
    den = den + pi(j+1);
end
if (M > m)
    % me parece que este if no es necesario
    for j=m:M-1
        den = den + pi(2*(j+1)+I);
    end
end


pi = StationaryMeasure(l,X1,lambda,mu,q,theta);

for j=0:X1
    for d=1:2
        num = num - pi(2*j+d) * cost2(j,b) ;
    end
end

m = min(l);
[M,I] = max(l);

for j=0:(2*m+1)
    den = den - pi(j+1);
end
if (M > m)
    % me parece que se puede poner el for sin el if
    for j=m:M-1
        den = den - pi(2*(j+1)+I);
    end
end

L = num / den ;

end

function [pi] = StationaryMeasure(n,X1,lambda,mu,q,theta)

%n = 2-dimension vector that indicates the threshold.
%X1 = where we truncate. In general, X1 = N + 50, where N is how many states
%we consider for a bandit.

if X1 <= max(n(1),n(2))
    disp('Truncation too low.')
    return
end
    

gamma = max(lambda(1),lambda(2)) + max(mu(1),mu(2)) + max(q(1),q(2)) + X1 * max(theta(1),theta(2)); % Uniformization parameter

pi = zeros(1,2*X1+2);
P = zeros(X1+1,2,X1+1,2); %anteultima columna = 1, ultima = 0.
P2 = zeros(2*X1+2,2*X1+2);
res = zeros(1,2*X1+2);
res(1,2*X1+1)=1;

for d=1:2
    for i=0:n(d)
        %up to n(d) there are no departures
        P(i+1,d,i+1,3-d) = q(d);
        P(i+1,d,i+2,d) = lambda(d);
        P(i+1,d,max(i,1),d) = i * theta(d);
        P(i+1,d,i+1,d) = gamma - q(d) - lambda(d) - i * theta(d);
    end
    for i=n(d)+1:X1
        P(i+1,d,i+1,3-d) = q(d);
        P(i+1,d,min(i+1,X1)+1,d) = lambda(d);
        P(i+1,d,max(i,1),d) = mu(d) + i * theta(d);
        P(i+1,d,i+1,d) = P(i+1,d,i+1,d) + gamma - q(d) - lambda(d) - mu(d) - i * theta(d);
    end
end


P = P/gamma;

for i=0:X1
    for d=1:2
        P(i+1,d,i+1,d) = P(i+1,d,i+1,d) - 1;
    end
end

for i =0:X1
    for d=1:2
        P(i+1,d,X1+1,1) = 1;
    end
end

for i=0:X1
    for d=1:2
        for j=0:X1
            for e=1:2
        P2(2*i+d,2*j+e)=P(i+1,d,j+1,e);
            end
        end
    end
end

pi = res / P2;

end

function f = cost2(j,b)
% We have linear cost

    f = j*b;

end
function sl = slope(n,X1,lambda,mu,q,theta)
pi = StationaryMeasure(n,X1,lambda,mu,q,theta);

sl = 0;
m = min(n);
[M,I] = max(n);

for j=0:(2*m+1)
    sl = sl + pi(j+1);
end
if (M > m)
    % me parece que este if no es necesario
    for j=m:M-1
        sl = sl + pi(2*(j+1)+I);
    end
end

sl = -sl;

end