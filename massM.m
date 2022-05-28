function [M] = massM(xh,shapeFn,noOfIntegPt)
    
    n = length(xh)-1;
    
    if shapeFn == 1
    
        M = sparse(n,n);   
        M(1,1) = meij(1,2,2,xh, shapeFn, noOfIntegPt);
    
        m1 = @(e)([meij(e,1,1,xh, shapeFn, noOfIntegPt) meij(e,1,2,xh, shapeFn, noOfIntegPt);
               0 meij(e,2,2,xh, shapeFn, noOfIntegPt)]); 
           
    for k = 2:n
        m = m1(k);
        m = m + tril(m',-1);
        M(k-1:k,k-1:k) = M(k-1:k,k-1:k) + m;
    end
           
    end
    
    if shapeFn == 2
        M = sparse(n*2,n*2); 
        %loads the inital 2x2
        M(1:2,1:2) = [meij(1,2,2,xh, shapeFn, noOfIntegPt) meij(1,2,3,xh, shapeFn, noOfIntegPt);
            0 meij(1,3,3,xh, shapeFn, noOfIntegPt)]; 
        M(2,1) = M(1,2);
        
        %a functional that creates a 3x3 Meij matrix that will be added to
        %the gloabl matrix
        m2 = @(e)([meij(e,1,1,xh, shapeFn, noOfIntegPt) meij(e,1,2,xh, shapeFn, noOfIntegPt) meij(e,1,3,xh, shapeFn, noOfIntegPt);
                    0 meij(e,2,2,xh, shapeFn, noOfIntegPt) meij(e,2,3,xh, shapeFn, noOfIntegPt);
                    0 0 meij(e,3,3,xh, shapeFn, noOfIntegPt)]);
       
        for k = 2:n
            %index
            s1 = 2*(k-1);
            s2 = 2*k;
            %calculates the uper diag matrix and then sets the uper part 
            %to the lower part
            m = m2(k);
            m = m + tril(m',-1);
            M(s1:s2,s1:s2) = M(s1:s2,s1:s2) + m;
        end

    end

end