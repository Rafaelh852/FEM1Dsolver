function [y] = fei(a, c, f, p0, e, i, xh, shapeFn, noOfIntegPt)

  x1 = xh(e);
  x2 = xh(e+1);  
  
  fn = @(x)((f(x)*shapeFn1d(i,x,x1,x2,shapeFn))...  
            -(a(x)*p0*shapeFnDer1d(1,x,xh(1),xh(2),shapeFn)*shapeFnDer1d(i,x,x1,x2,shapeFn))... 
            -(c(x)*p0*shapeFn1d(1,x,xh(1),xh(2),shapeFn)*shapeFn1d(i,x,x1,x2,shapeFn)));

  y = gaussQuad1d(fn,x1,x2, noOfIntegPt);
end