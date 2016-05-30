%% 2-D Case
clear
clc

% Create symbolic state vector
x = sym('x',[2,1]) ;

% Create paraboloid function f, and compute the gradient analytically
fsym = 1 - x(1)^2 - x(2)^2 ;
fMultiIn = matlabFunction(fsym) ;
gMultiIn = matlabFunction(gradient(fsym)) ;

% Convert the symbolic f and g to anonymous functions
f = @(x) callFunction(fMultiIn,x) ;
g = @(x) callFunction(gMultiIn,x) ;
    % Note that 'callFunction.m' contains the following:
    % function out = callFunction(f,x)
    %     x = num2cell(x) ;
    %     out = f(x{:}) ;
    % end

% Create matrix of points to test, where each column is a point
X = [1 0 -1 0;
     0 1 0 -1];

% Compare numerical and analytical gradients
gradientComparison(f,g,X)

%% n-D Generic Case
clear
clc

n = 3 ; % number of states
N = 5 ; % number of points to evaluate
k = 5 ; % max power of any individual state; k*n is the max order possible

% Create symbolic state vector
x = sym('x',[n,1]) ;

% Generate polynomial symbolic function f, and compute gradient analytically
fsym = 0 ;
for idx1 = 1:n
    temp = 1 ;
    for idx2 = 1:n
        temp = temp.*x(idx2).^(randi(k,1)) ;
    end
    fsym = fsym + temp ;
end
gsym = gradient(fsym) ;

fMultiIn = matlabFunction(fsym) ;
gMultiIn = matlabFunction(gsym) ;

% Convert the symbolic f and g to anonymous functions
f = @(x) callFunction(fMultiIn,x) ;
g = @(x) callFunction(gMultiIn,x) ;

% Create matrix of points to test, where each column is a point
X = 2.*rand(n,N) ;

% Compare numerical and analytical gradients
[err,okFlag] = gradientComparison(f,g,X) ;
disp(err)

%% Vector Fixed Case
clear
clc

syms x1 x2 x3 ;

fsym = [-x1^2+x2^2; x3^3 - x2^2] ;
gsym = jacobian(fsym) ;
fMultiIn = matlabFunction(fsym) ;
gMultiIn = matlabFunction(gsym) ;
f = @(x) callFunction(fMultiIn,x) ;
g = @(x) callFunction(gMultiIn,x) ;

v = [1 1 1]' ;

gnum = numjacobian(f,v) ;
ganl = g(v) ;
err = norm(gnum - ganl).^2 ;
disp(err)

%% Vector Random Case
clear
clc

n = 2 ; % state dimension
m = 2 ; % vector function dimension
k = 3 ; % max power of any individual state; k*n is the max order possible

% Create symbolic state vector
x = sym('x',[n,1]) ;

% Generate symbolic function f, which has negative exponential dynamics;
% the numerical and analytic jacobians tidily converge to one another
fsym = sym(ones(m,1)) ;
for idx1 = 1:m
    for idx2 = 1:n
        fsym(idx1) = fsym(idx1).*exp(-x(idx2).^(randi(k,1))) ;
    end
end
gsym = jacobian(fsym) ;

fMultiIn = matlabFunction(fsym) ;
gMultiIn = matlabFunction(gsym) ;
f = @(x) callFunction(fMultiIn,x) ;
g = @(x) callFunction(gMultiIn,x) ;

err = 0 ;
for test = 1:100
v = 1.*rand(n,1) ;

gnum = numjacobian(f,v) ;
ganl = g(v) ;
err = max(err,norm(gnum - ganl).^2) ;
end

disp(err)

%% testing the better numerical jacobian
clear
clc

fxup = @(x,u,p) [x(1)*u(1)*p(1); x(2)*u(1)*p(1)] ;
gxp = @(x,u,p) [u(1)*p(1), 0, x(1)*u(1);
                0, u(1)*p(1), x(2)*u(1)];
in = {rand(2,1),rand(1),rand(1)} ;
J = numericJacobian(fxup,[1 3],in{:}) ;
gxp(in{:}) ;

%% testing numerical Jacobian with test functions
clear
clc

in = {rand(4,1),rand(1),rand(3,1)} ;
[dxdt, dfdx, dfdp] = testNumJacob(in{:}) ;
Janl = [dfdx, dfdp] ;
Jnum = numericJacobian(@testNumJacob,[1 3],in{:}) ;
norm(Janl - Jnum)


