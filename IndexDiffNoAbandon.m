function[I]=IndexDiffNoAbandon(m,mu)

I = zeros(length(m),length(m),2,2); % Initializing the array that keeps Whittle Index
TrunI = zeros(length(m)-1,length(m)-1,2,2);
Iddle1 = zeros(length(m)-1,1) + 1;
Iddle2 = zeros(1,length(m)-1) + 2;
for d=1:2
    for e=1:2
        if mu(d,1) > mu(3-d,1) && mu(e,2) < mu(3-e,2)
        %IF 1 IS IN GOOD ENVIR AND 2 IS IN BAD ENVIR
            TrunI(:,:,d,e) = zeros(length(m)-1,length(m)-1) + 1;
        elseif mu(d,1) < mu(3-d,1) && mu(e,2) > mu(3-e,2)
        %IF 1 IS IN BAD ENVIR AND 2 IS IN GOOD ENVIR
            TrunI(:,:,d,e) = zeros(length(m)-1,length(m)-1) + 2;
        else
            for j=1:(length(m)-1)
                for k=1:(length(m)-1)
                    if (m(1,j+1,d)<=m(2,k+1,e)) 
                    TrunI(j,k,d,e) = 2;
                    elseif (m(1,j+1,d)>m(2,k+1,e))
                    TrunI(j,k,d,e) = 1;
                    end
                end
            end
        end
       I(:,:,d,e) = [1 Iddle2 ; Iddle1 TrunI(:,:,d,e)]; 
    end
end

end
