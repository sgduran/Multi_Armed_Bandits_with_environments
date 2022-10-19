%CODE FOR GETTING INDICES GIVEN THE PARAMETERS

lambda = [4 4; 4 4];
mu = [8 27; 5 21];
theta = [0.1 0.4; 0.1 0.3];
%q = [15 17];
qD = [15 17; 10 2];

%W = WIndices(mu,theta,q);

WD = WIndicesD(mu,theta,qD);

%W;
WD


function [W] = WIndices(mu,theta,q)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+q(1)+q(2)) / (theta(1,i)*theta(2,i)+theta(2,i)*q(1)+theta(1,i)*q(2)));
    end
end

end

function [W] = WIndicesD(mu,theta,qD)

W = zeros(2,2);

for i=1:2
    for p=1:2
        W(p,i) = (mu(p,i)) * ((theta(3-p,i)+qD(i,1)+qD(i,2)) / (theta(1,i)*theta(2,i)+theta(2,i)*qD(i,1)+theta(1,i)*qD(i,2)));
    end
end

end