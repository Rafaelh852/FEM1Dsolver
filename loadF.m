function [F] = loadF(a, c, f, p0, QL, xh, shapeFn, noOfIntegPt)

    n = length(xh)-1;
    
    
    if shapeFn == 1
        
        F = zeros(n,1);
        F(1) = fei(a, c, f, p0, 1, 2, xh, shapeFn, noOfIntegPt) + QL;
        
        f1 = @(e)([fei(a, c, f, p0, e, 1, xh, shapeFn, noOfIntegPt)-QL;
                   fei(a, c, f, p0, e, 2, xh, shapeFn, noOfIntegPt)+QL]);
        
        for k = 2:n
           F(k-1:k) = F(k-1:k) + f1(k);  
        end
        
    end
    
    if shapeFn == 2
        F = zeros(n*2,1);
        
        F(1:2) = [fei(a, c, f, p0, 1, 2, xh, shapeFn, noOfIntegPt);
                  fei(a, c, f, p0, 1, 3, xh, shapeFn, noOfIntegPt)+QL];
          
        f2 = @(e)([fei(a, c, f, p0, e, 1, xh, shapeFn, noOfIntegPt)-QL;
                   fei(a, c, f, p0, e, 2, xh, shapeFn, noOfIntegPt);
                   fei(a, c, f, p0, e, 3, xh, shapeFn, noOfIntegPt)+QL]);
              
        for k = 2:n
            s1 = 2*(k-1);
            s2 = 2*k;
            F(s1:s2) = F(s1:s2) + f2(k); 
        end
        
    end

end