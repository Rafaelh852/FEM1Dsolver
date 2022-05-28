% these are the linear and quadratic shape functions on a 1-d element [x1,x2]
% y = shapeFn1d(i,x,x1,x2,p)
%
% input:
%   i: i-th shape function
%   x: input for shape function
%   x1, x2: end points for the 1d element, x1 < x2
%   p: order of shape function
% output:
%   y: thise shapeFn(x)
%
% Last update: Chung-min Lee Mar 4, 2017

function y = shapeFn1d(i,x,x1,x2,p)
ind = (x < x2+1e-12 & x > x1-1e-12);
y = zeros(size(x));

%linear shape functions
if (p == 1)
    if (i == 1)
        y(ind) = (x2(ind) - x(ind))./(x2(ind)-x1(ind));
    elseif (i == 2)
        y(ind) = (x(ind) - x1(ind))./(x2(ind) - x1(ind));
    else
        y = [];
    end
end

% quadratic shape functions
if (p == 2)
    
    if ( i == 1)
        y(ind) = (x2(ind) - x(ind)).*(x2(ind)-2*x(ind)+x1(ind))./(x2(ind)-x1(ind)).^2;
    elseif ( i == 2)
        y(ind) = 4*(x(ind)-x1(ind)).*(x2(ind)-x(ind))./(x2(ind)-x1(ind)).^2;
    elseif ( i ==3)
        y(ind) = (x1(ind)-x(ind)).*(x2(ind)-2*x(ind)+x1(ind))./(x2(ind)-x1(ind)).^2;
    else
        y = [];
    end
end

if (p ~= 1 && p ~= 2)
    y = [];
end