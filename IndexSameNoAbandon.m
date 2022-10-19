function[I]=IndexSameNoAbandon(m,mu)

I = zeros(length(m),length(m),2); % Initializing the array that keeps Whittle Index
TrunI = zeros(length(m)-1,length(m)-1,2);
Iddle1 = zeros(length(m)-1,1) + 1;
Iddle2 = zeros(1,length(m)-1) + 2;
for d=1:2
    if mu(d,1) > mu(3-d,1) && mu(d,2) < mu(3-d,2)
    %IF 1 IS IN GOOD ENVIR AND 2 IS IN BAD ENVIR
        TrunI(:,:,d) = zeros(length(m)-1,length(m)-1) + 1;
    elseif mu(d,1) < mu(3-d,1) && mu(d,2) > mu(3-d,2)
    %IF 1 IS IN BAD ENVIR AND 2 IS IN GOOD ENVIR
        TrunI(:,:,d) = zeros(length(m)-1,length(m)-1) + 2;
    else
        for j=1:(length(m)-1)
            for k=1:(length(m)-1)

                if (m(1,j+1,d)<=m(2,k+1,d)) 
                TrunI(j,k,d) = 2;
                elseif (m(1,j+1,d)>m(2,k+1,d))
                TrunI(j,k,d) = 1;
                end
            end

        end
    end
    I(:,:,d) = [1 Iddle2 ; Iddle1 TrunI(:,:,d)];
end



end
