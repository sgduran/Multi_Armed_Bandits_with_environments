function[I]=IndexMatrixDiff(m)

I = zeros(length(m),length(m),2,2); % Initializing the array that keeps Whittle Index
for j=0:(length(m)-1)
    for k=0:(length(m)-1)
        for d=1:2
            for e=1:2
                if (m(1,j+1,d)<=m(2,k+1,e)) 
                I(j+1,k+1,d,e) = 2;
                elseif (m(1,j+1,d)>m(2,k+1,e))
                I(j+1,k+1,d,e) = 1;
                end
            end
        end
        
    end
end

end
