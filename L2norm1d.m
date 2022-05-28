function [y] = L2norm1d(f,a,b,noOfEle)
    y = 0;
    fn = @(x)(f(x)*f(x)); 
    step = (b-a)/noOfEle;
    xh = [a:step:b];
    
    for k = 1:noOfEle
      x1 = xh(k);
      x2 = xh(k+1);
      y = y + gaussQuad1d(fn,x1,x2,3);
    end
    
    y = sqrt(y);
end
