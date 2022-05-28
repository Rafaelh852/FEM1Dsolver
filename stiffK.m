function [K] = stiffK(a, c, xh, shapeFn, noOfIntegPt)

    n = length(xh)-1;
    
    if shapeFn == 1
    
    K = sparse(n,n);
    K(1,1) = keij(a,c,1,2,2,xh, shapeFn, noOfIntegPt);
    
    k1 = @(e)([keij(a,c,e,1,1,xh, shapeFn, noOfIntegPt) keij(a,c,e,1,2,xh, shapeFn, noOfIntegPt);
               0 keij(a,c,e,2,2,xh, shapeFn, noOfIntegPt)]); 
           
    for k = 2:n
        m = k1(k);
        m = m + tril(m',-1);
        K(k-1:k,k-1:k) = K(k-1:k,k-1:k) + m;
    end
           
    end
    
    if shapeFn == 2
        
        K = sparse(n*2,n*2);
        
        %loads the inital 2x2
        K(1:2,1:2) = [keij(a,c,1,2,2,xh, shapeFn, noOfIntegPt) keij(a,c,1,2,3,xh, shapeFn, noOfIntegPt);
            0 keij(a,c,1,3,3,xh, shapeFn, noOfIntegPt)]; 
        K(2,1) = K(1,2);
        
        %a functional that creates a 3x3 Meij matrix that will be added to
        %the gloabl matrix
        k2 = @(e)([keij(a,c,e,1,1,xh, shapeFn, noOfIntegPt) keij(a,c,e,1,2,xh, shapeFn, noOfIntegPt) keij(a,c,e,1,3,xh, shapeFn, noOfIntegPt);
                    0 keij(a,c,e,2,2,xh, shapeFn, noOfIntegPt) keij(a,c,e,2,3,xh, shapeFn, noOfIntegPt);
                    0 0 keij(a,c,e,3,3,xh, shapeFn, noOfIntegPt)]);
       
        for k = 2:n
            %index
            s1 = 2*(k-1);
            s2 = 2*k;
            %calculates the uper diag matrix and then sets the uper part 
            %to the lower part
            m = k2(k);
            m = m + tril(m',-1);
            K(s1:s2,s1:s2) = K(s1:s2,s1:s2) + m;
        end

    end


end