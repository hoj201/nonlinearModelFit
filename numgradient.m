function g = numgradient(f,x)
% For a function f, compute the numerical gradient at point or points x,
% where x is either a column vector or matrix of concatenated column
% vectors, with a number of rows equivalent to the dimensionality of x

[n,N] = size(x) ; % n is dimension of x, N is number of points to evaluate

h = max(sqrt(eps),sqrt(eps)*norm(x)) ; % step size
hmat = h.*eye(n) ;
g = zeros(n,N) ; % gradient matrix to fill in

for idxCol = 1:N
    for idxRow = 1:n
        f0 = f(x(:,idxCol) - hmat(:,idxRow)) ;
        f1 = f(x(:,idxCol) + hmat(:,idxRow)) ;
        g(idxRow,idxCol) = (f1 - f0)./(2*h) ;
    end
end

end
