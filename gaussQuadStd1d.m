function [y] = gaussQuadStd1d(g,noOfIntegPt)
    
    if noOfIntegPt ==2 
        
        y = g(-sqrt(3)/3) + g(sqrt(3)/3);
        
       return 
    end

    if noOfIntegPt == 3
        y = (5/9)*g(-sqrt(3/5)) + (8/9)*g(0) + (5/9)*g(sqrt(3/5));
        
        return
    end
    
    "noOfIntegPt is not 2 or 3"
    y = 0;
    return
end
