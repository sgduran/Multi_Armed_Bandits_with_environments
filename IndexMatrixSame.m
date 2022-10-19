function[I]=IndexMatrixSame(m)

I = zeros(length(m),length(m),2); % Initializing the array that keeps Whittle Index
for j=0:(length(m)-1)
    for k=0:(length(m)-1)
        for d=1:2
            if (m(1,j+1,d)<=m(2,k+1,d)) 
            I(j+1,k+1,d) = 2;
            elseif (m(1,j+1,d)>m(2,k+1,d))
            I(j+1,k+1,d) = 1;
            end
        end
        
    end
end

end
