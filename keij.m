function [y] =  keij(a, c, e, i, j, xh, shapeFn, noOfIntegPt)
  x1 = xh(e);
  x2 = xh(e+1);
  
  fn = @(x)(a(x)*shapeFnDer1d(i,x,x1,x2,shapeFn)*shapeFnDer1d(j,x,x1,x2,shapeFn) + c(x)*shapeFn1d(i,x,x1,x2,shapeFn)*shapeFn1d(j,x,x1,x2,shapeFn));
  
  y = gaussQuad1d(fn,x1,x2, noOfIntegPt);
  
end