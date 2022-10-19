function[I]=IndexMatrixRandomSame(N1)

I = zeros(N1+1,N1+1,2); % Initializing the array that keeps Whittle Index
I(:,1,:) = ones(N1+1,2);
I(1,:,:) = 2*ones(N1+1,2);
I(2:end,2:end,:) = randi(2,N1,N1,2);

end
