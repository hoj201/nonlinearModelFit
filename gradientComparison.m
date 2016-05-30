function [err, okFlag] = gradientComparison(f,g,x)
% Given a function f and an analytic gradient g, compute the numerical
% gradient of f about point or points x and compare it with g. Return the
% error between the numerical and and analytic gradients. The parameter x
% is assumed to be either a column vector or a matrix of concatenated
% column vectors, with a number of rows equal to the dimension of x.

[n,N] = size(x) ; % n is dimension of x, N is number of points to evaluate

% compute numerical gradient
gnum = numgradient(f,x) ;

% compute analytic gradient at each point
ganl = zeros(n,N) ; % gradient matrix to fill in
for idx = 1:N
    ganl(:,idx) = g(x(:,idx)) ;
end

err = max(sum((gnum - ganl).^2,2)) ;
okFlag = err < 1e-6 ;

end