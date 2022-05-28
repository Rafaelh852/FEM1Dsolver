function [uh]  = myFE1dibvp(a, c, f, p0, QL, u0, L, T, dt, noOfEle, shapeFn)
    

    s = noOfEle*shapeFn;
    dx = L/noOfEle;
    xh = [0:dx:L];
    th = [0:dt:T];
    
    p = shapeFn;
    q = shapeFn;
    
    M = massM(xh,shapeFn,3);
    K = stiffK(a,c, xh, shapeFn, 3);

    f1 = @(t)(@(x)(f(x,t)));
    F = @(t)(loadF(a, c, f1(t), p0, QL(t), xh, shapeFn, 3));
    
    if shapeFn == 2
        dx = L/s;
        xh = [0:dx:L];
    end
  
    W = zeros(s,length(th));
    
    LHS = (M +.5*dt.*K);
    LU = bandLU(LHS,p,q) ;
    lower = tril(LU,-1) + eye(length(LU));
    upper = (tril(LU'))';
    
    for t = 1:length(th)
        if t == 1
         for e = 2:s+1
           W(e-1,t) = u0(xh(e));
         end
        else
           t1 = th(t-1);
           t2 = th(t);
           RHS1 = ((M-.5*dt*K)*W(:,t-1))+.5*dt*(F(t1) + F(t2));
           RHS2 = bandforward(lower,RHS1,p);
           W(:,t) = bandbackward(upper,RHS2,p);
            
        end
        
        uh{t} = approxSol(W(:,t), p0, xh(1:shapeFn:end), shapeFn); 
    end
    
    %debugging purposes uncomment W to check vals
    %W
    
end

function [K] = bandLU(K,p,q)

% test K
%K= [1 2 0 0;3 1 2 2;0 0 1 3;0 0 1 1];
%[p,q] = bandwidth(K);

n = length(K);

for m = 1:n-1
    for i = m+1:min([m+p,n])
        K(i,m) = K(i,m)/K(m,m);
    end
    for j = m+1:min([m+q,n])
        for i = m+1:min([m+p,n])
            K(i,j)= K(i,j) - K(i,m)*K(m,j);
        end
    end
end

K = K;
% %original K == (tril(K,-1) + eye(n))*triu(K)
% K
% L = tril(K,-1) + eye(n)
% U = triu(K)

% %test vector
% f = [2; 2; 3; 4];
% 
% %--- Bandforward ---
% for j = 1:n
%     for i= j+1:min([j+p,n])
%         f(i) = f(i) - L(i,j)*f(j);
%     end 
% end
% 
% %checking solution
% f
% L*f
% 
% %--- Bandbackward ---    
% for j = n:-1:1
%     f(j) = f(j)/U(j,j);
%     for i  = max([1,j-q]):j-1
%         f(i) = f(i) - U(i,j)*f(j);
%     end
% end
% 
% %checking solution
% f
% U*f

end

function [z] = bandbackward(L,f,p)

n = length(f);
  
for j = n:-1:1
    f(j) = f(j)/L(j,j);
    for i  = max([1,j-p]):j-1
        f(i) = f(i) - L(i,j)*f(j);
    end
end
z = f;
end

function [z] = bandforward(L,f,p)

n = length(f);

for j = 1:n
    for i= j+1:min([j+p,n])
        f(i) = f(i) - L(i,j)*f(j);
    end 
end
 z = f;
end

function [y] = shapeFn1d(i,x,x1,x2,p)

    if x <x1 || x > x2 
      % "x is not in [x1,x2]"
      y = 0;
      return 
    end

    h = x2-x1;
    % linear shape function
    if p==1
        if i == 1
            y = (x2 - x)/h;
            return
        end
        
        if i == 2 
            y = (x - x1)/h;
            return
        end
        "i  is not 1 or 2"
        y=0;
      return
    end
    
    % quadratic shape function
    if p == 2
        
        if i == 1
            y = (1/(h*h))*(x2-x)*(x2-2*x+x1);
            return
        end
        
        if i == 2
            y = (4/(h*h))*(x-x1)*(x2-x);
            return
        end
         
        if i == 3
            y = (1/(h*h))*(x1-x)*(x2-2*x+x1);
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

function [y] = gaussQuad1d(fn,lowerLimit,upperLimit, noOfIntegPt)
 b = upperLimit; 
 a = lowerLimit;
 g = @(x)(fn(b +(b-a)*(x-1)/2));
 y = ((b-a)/2)*gaussQuadStd1d(g,noOfIntegPt);

end

function [y] = meij(e,i,j,xh, shapeFn, noOfIntegPt)
 
  x1 = xh(e);
  x2 = xh(e+1);
  
  fn = @(x)(shapeFn1d(i,x,x1,x2,shapeFn)*shapeFn1d(j,x,x1,x2,shapeFn));
 
  y = gaussQuad1d(fn,x1,x2, noOfIntegPt);
  

end

function [y] =  keij(a, c, e, i, j, xh, shapeFn, noOfIntegPt)
  x1 = xh(e);
  x2 = xh(e+1);
  
  fn = @(x)(a(x)*shapeFnDer1d(i,x,x1,x2,shapeFn)*shapeFnDer1d(j,x,x1,x2,shapeFn) + c(x)*shapeFn1d(i,x,x1,x2,shapeFn)*shapeFn1d(j,x,x1,x2,shapeFn));
  
  y = gaussQuad1d(fn,x1,x2, noOfIntegPt);
  
end

function [y] = fei(a, c, f, p0, e, i, xh, shapeFn, noOfIntegPt)

  x1 = xh(e);
  x2 = xh(e+1);  
  
  fn = @(x)((f(x)*shapeFn1d(i,x,x1,x2,shapeFn))...  
            -(a(x)*p0*shapeFnDer1d(1,x,xh(1),xh(2),shapeFn)*shapeFnDer1d(i,x,x1,x2,shapeFn))... 
            -(c(x)*p0*shapeFn1d(1,x,xh(1),xh(2),shapeFn)*shapeFn1d(i,x,x1,x2,shapeFn)));

  y = gaussQuad1d(fn,x1,x2, noOfIntegPt);
  %y = gaussQuad1d(fn,x1,x2, noOfIntegPt) + Q(x2)*shapeFn1d(i,x2,x1,x2,shapeFn)- Q(x1)*shapeFn1d(i,x1,x1,x2,shapeFn);
end

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
        M = sparse(2*n,2*n); 
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

function [F] = loadF(a, c, f, p0, QL, xh, shapeFn, noOfIntegPt)

    n = length(xh)-1;
    
    
    if shapeFn == 1
        
        F = zeros(n,1);
        F(1) = fei(a, c, f, p0, 1, 2, xh, shapeFn, noOfIntegPt);
        F(n) = QL;
        
        f1 = @(e)([fei(a, c, f, p0, e, 1, xh, shapeFn, noOfIntegPt);
                   fei(a, c, f, p0, e, 2, xh, shapeFn, noOfIntegPt)]);
        
        for k = 2:n
           F(k-1:k) = F(k-1:k) + f1(k);  
        end
        
    end
    
    if shapeFn == 2
        F = zeros(n*2,1);
        
        F(1:2) = [fei(a, c, f, p0, 1, 2, xh, shapeFn, noOfIntegPt);
                  fei(a, c, f, p0, 1, 3, xh, shapeFn, noOfIntegPt)];
        F(n*2) = QL;     
        f2 = @(e)([fei(a, c, f, p0, e, 1, xh, shapeFn, noOfIntegPt);
                   fei(a, c, f, p0, e, 2, xh, shapeFn, noOfIntegPt);
                   fei(a, c, f, p0, e, 3, xh, shapeFn, noOfIntegPt)]);
              
        for k = 2:n
            s1 = 2*(k-1);
            s2 = 2*k;
            F(s1:s2) = F(s1:s2) + f2(k); 
        end
        
    end

end

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

