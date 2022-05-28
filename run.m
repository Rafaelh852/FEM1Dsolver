clear
a = @(x) 6-2*x;
c = @(x) ((2*pi*pi)*x);
f = @(x,t) ((-x-2*t+6*pi*pi*(4-x*t))*cos(pi*x)-2*pi*(6*t-3*x*t+4)*sin(pi*x));
p0= 4;
QL = @(t)(-2*t);
u0 = @(x)(4*cos(pi*x));
L=2;
T=1;

% a = @(x) 4-x;
% c = @(x) x;
% f = @(x,t) -2*t.*sin(x)/((1+t.^2).^2) + (cos(x)+4*sin(x))./(1+t.^2);
% L = 2;
% T = 1;
% p0 = 0;
% QL = @(t) 2*cos(2)/(1+t.^2);
% u0 =@(x) sin(x);

dt = 0.01;
% shapeFn = 2;
% noOfEle = 4;

%exact solu
u = @(x,t)((4-x.*t).*cos(pi*x));

%aprox solu
% uh  = myFE1dibvp(a, c, f, p0, QL, u0, L, T, dt, noOfEle, shapeFn);
uht  =@(e,s)( myFE1dibvp(a, c, f, p0, QL, u0, L, T, dt, e, s));

uh16linear = uht(16,1);
uh16quad = uht(16,2);

uh32linear = uht(32,1);
uh32quad = uht(32,2);

uh64linear = uht(64,1);
uh64quad = uht(64,2);

uh128linear = uht(128,1);
uh128quad = uht(128,2);

[x16p y16p z16p] = plotuh(uh16linear,L,T,dt,16,p0);
[x16q y16q z16q] = plotuh(uh16quad,L,T,dt,16,p0);

[x32p y32p z32p] = plotuh(uh32linear,L,T,dt,32,p0);
[x32q y32q z32q] = plotuh(uh32quad,L,T,dt,32,p0);

[x64p y64p z64p] = plotuh(uh64linear,L,T,dt,64,p0);
[x64q y64q z64q] = plotuh(uh64quad,L,T,dt,64,p0);

[x128p y128p z128p] = plotuh(uh128linear,L,T,dt,128,p0);
[x128q y128q z128q] = plotuh(uh128quad,L,T,dt,128,p0);

zu =@(X,Y)((4-X.*Y).*cos(pi*X));

z128e = zu(x128p,y128p);

figure
subplot(3,1,1);
plot(x16p(:,end),z16p(:,end),x32p(:,end),z32p(:,end),x64p(:,end),z64p(:,end),x128p(:,end),z128p(:,end),x128p(:,end),z128e(:,end));
title("approximation using linear shape functions");
legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128","exact solution");
xlabel("x");
ylabel("approximated sol.");

subplot(3,1,2);
plot(x16q(:,end),z16q(:,end),x32q(:,end),z32q(:,end),x64q(:,end),z64q(:,end),x128q(:,end),z128q(:,end),x128p(:,end),z128e(:,end));
title("approximation using quadratic shape functions");
legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128","exact solution");
xlabel("x");
ylabel("approximated sol.");

subplot(3,1,3);
plot(x128p(:,end),z128p(:,end),x128q(:,end),z128q(:,end),x128p(:,end),z128e(:,end));
title("aproximation with 16 elements");
legend("linear shape function","quadratic shape function","exact solution");
xlabel("x");
ylabel("approximated sol.");

z16e = zu(x16p,y16p);
z32e = zu(x32p,y32p);
z64e = zu(x64p,y64p);

err1p = err(z16p,z16e,0);
err1q = err(z16q,z16e,0);
err2p = err(z32p,z32e,0);
err2q = err(z32q,z32e,0);
err3p = err(z64p,z64e,0);
err3q = err(z64q,z64e,0);
err4p = err(z128p,z128e,0);
err4q = err(z128q,z128e,0);

figure
subplot(3,1,1);
semilogy(x16p(:,end),err1p(:,end),x32p(:,end),err2p(:,end),x64p(:,end),err3p(:,end),x128p(:,end),err4p(:,end));
title("abs. error using linear shape functions");
legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128");
xlabel("x");
ylabel("abs error");

subplot(3,1,2);
semilogy(x16q(:,end),err1q(:,end),x32q(:,end),err2q(:,end),x64q(:,end),err3q(:,end),x128q(:,end),err4q(:,end));
title("abs. error using quadratic shape functions");
legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128");
xlabel("x");
ylabel("abs error");

subplot(3,1,3);
semilogy(x128p(:,end),err4p(:,end),x128q(:,end),err4q(:,end));
title("abs. error using 128 elements");
legend("linear shape function","quadratic shape function");
xlabel("x");
ylabel("abs error");

figure
subplot(2,1,1)
mesh(x128p,y128p,err4p);
title("abs. error between true sol and linear shape fn with 128 elements");
xlabel("x");
ylabel("t");
zlabel("abs error");

subplot(2,1,2)
mesh(x128p,y128p,err4q);
title("abs. error between true sol and linear shape fn with 128 elements");
xlabel("x");
ylabel("t");
zlabel("abs error");


% finding convergence
shapeFn = 1;
noOfEle = 16;
LM(:,1) = ones(3,1);
LM(1,2) = log(1/8);
LM(2,2) = log(1/16);
LM(3,2) = log(1/32);
utrue = @(x)((4-x*T)*cos(pi*x));

for i = 1:3
    
    uh1 = uht(noOfEle,shapeFn);
    fL2 = @(x)(abs(utrue(x)-uh1{end}(x)));
    
    L2 = L2norm1d(fL2,0,L,noOfEle);
    RV(i,1) = log(L2);
    
    noOfEle = 2*noOfEle;
end

CV1 = LM\RV;
c1 = CV1(2);

shapeFn = 2;
noOfEle = 16;
LM2(:,1) = ones(4,1);
LM2(1,2) = log(1/8);
LM2(2,2) = log(1/16);
LM2(3,2) = log(1/32);
LM2(4,2) = log(1/64);

for i = 1:4
    
    uh2 = uht(noOfEle,shapeFn);
    fL2 = @(x)(abs(utrue(x)-uh2{end}(x)));
    
    L2 = L2norm1d(fL2,0,L,noOfEle);
    RV2(i,1) = log(L2);
    
    noOfEle = 2*noOfEle;
end

CV2 = LM2\RV2;
c2 = CV2(2);


function [X Y Z] = plotuh(uh,L,T,dt,noOfEle,p0)
%     dx = L/(2*noOfEle);
%     x = 0:dx:L;
    x = 0:0.01:L;
    t = 0:dt:T;
    n = length(x);
    m = length(t);
    
    [X Y] = meshgrid(x,t);
     X = X';
     Y = Y';
    
    Z = zeros(n,m);
    
    for j = 1:m
        Z(1,j) = p0;
        for i = 2:n
            Z(i,j) = uh{j}(x(i));
        end
    end
   
end

function [ernum] = err(A,B,p)
    diff = abs(B-A);
    if p  == 0
        ernum = diff;
        
    else
       ernum = norm(abs(diff),p);
        
    end
        
end

%at time t where 0 < t < T 
% t = T;
% times = 0:dt:T;
% ts = find(times==t);
% 
% x = 0:.01:L;
% t = 0:dt:T;
% 
% figure
% subplot(3,1,1);
% hold on;
% plot(x,uh16linear{end}(x));
% plot(x,uh32linear{end}(x));
% plot(x,uh64linear{end}(x));
% plot(x,uh128linear{end}(x));
% plot(x,u(x,T));
%  title("approximation using linear shape functions");
% legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128", "exact solution");
% xlabel("x");
% ylabel("approximated sol.");
% 
% subplot(3,1,2);
% hold on;
% plot(x,uh16quad{end}(x));
% plot(x,uh32quad{end}(x));
% plot(x,uh64quad{end}(x));
% plot(x,uh128quad{end}(x));
% plot(x,u(x,T));
% title("approximation using quadratic shape functions");
% legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128","exact solution");
% xlabel("x");
% ylabel("approximated sol.");
% 
% subplot(3,1,3);
% hold on;
% plot(x,uh128linear{end}(x));
% plot(x,uh128quad{end}(x));
% plot(x,u(x,T));
% title("aproximation with 128 elements");
% legend("linear shape function","quadratic shape function","exact solution");
% xlabel("x");
% ylabel("approximated sol.");
% 
% 
% 
% figure
% subplot(3,1,1);
% hold on
% semilogy(x,abs(u(x,T) - uh16linear{end}(x)));
% semilogy(x,abs(u(x,T) - uh32linear{end}(x)));
% semilogy(x,abs(u(x,T) - uh64linear{end}(x)));
% semilogy(x,abs(u(x,T) - uh128linear{end}(x)));
% title("abs. error using linear shape functions");
% legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128");
% xlabel("x");
% ylabel("abs error");
% 
% 
% subplot(3,1,2);
% hold on
% semilogy(x,abs(u(x,T) - uh16quad{end}(x)));
% semilogy(x,abs(u(x,T) - uh32quad{end}(x)));
% semilogy(x,abs(u(x,T) - uh64quad{end}(x)));
% semilogy(x,abs(u(x,T) - uh128quad{end}(x)));
% title("abs. error using quadratic shape functions");
% legend("noOfEle = 16","noOfEle = 32","noOfEle = 64","noOfEle = 128");
% xlabel("x");
% ylabel("abs error");

% subplot(3,1,3);
% hold on
% semilogy(x,abs(u(x,T) - uh128linear{end}(x)));
% semilogy(x,abs(u(x,T) - uh128quad{end}(x)));
% title("abs. error using 128 elements");
% legend("linear shape function","quadratic shape function");
% xlabel("x");
% ylabel("abs error");
% 
% [X Y] = meshgrid(x,t);
%  
% for ts = 1:length(t)
%     Zlinear(ts,:) = uh128linear{ts}(x);
%     Zquad(ts,:) = uh128quad{ts}(x);
% end
% 
% figure
% subplot(1,2,1)
% mesh(X',Y',abs(u(X',Y') - Zlinear'));
% title("abs. error between true sol and linear shape fn with 128 elements");
% xlabel("x");
% ylabel("t");
% zlabel("abs error");
% 
% subplot(1,2,2)
% mesh(X',Y',abs(u(X',Y') - Zquad'));
% title("abs. error between true sol and linear shape fn with 128 elements");
% xlabel("x");
% ylabel("t");
% zlabel("abs error");
