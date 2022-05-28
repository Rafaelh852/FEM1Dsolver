function [y] = shapeFnDer1d(i,x,x1,x2,p)

    if x <x1 || x > x2 
      % "x is not in [x1,x2]"
      y = 0;
      return 
    end

    h = x2-x1;
    % linear shape function
    if p==1
        if i == 1
            y = -1/h;
            return
        end
        
        if i == 2 
            y = 1/h;
            return
        end
        "i  is not 1 or 2"
        y=0;
      return
    end
    
    % quadratic shape function
    if p == 2
        
        if i == 1
            y = (1/(h*h))*((x2-x)*(-2)-(x2-2*x+x1));
            return
        end
        
        if i == 2
            y = (4/(h*h))*((x2-x)-(x-x1));
            return
        end
         
        if i == 3
            y = (1/(h*h))*((x1-x)*(-2)-(x2-2*x+x1));
            return
        end
        
        "i is not 1,2, or 3"
        y = 0;
     return
    end
    
    "p is not 1 or 2"
    y=0;
    return



end