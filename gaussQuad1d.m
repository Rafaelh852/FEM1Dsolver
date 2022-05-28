function [y] = gaussQuad1d(fn,lowerLimit,upperLimit, noOfIntegPt)
 b = upperLimit; 
 a = lowerLimit;
 g = @(x)(fn(b +(b-a)*(x-1)/2));
 y = ((b-a)/2)*gaussQuadStd1d(g,noOfIntegPt);

end
