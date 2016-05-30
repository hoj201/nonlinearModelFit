function J = numjacobian(f,x)
% For a vector function f, of size Nf x 1, and a datum x, compute the
% numerical Jacobian matrix such that
%
%           J(i,j) = partial(f_i)/partial(X_j)
%
% Note that, in general, Nf = Nx, but this code will work for a vector
% function f that is of different dimension than its domain
%
% INPUTS:
%   f  -  Ndynmcs x 1 vector function
%   x  -  Nstates x 1 numeric vector
%
% OUTPUTS:
%   J  -  matrix of size Nf x Nx

Nx = length(x) ; % determine number of states
Nf = length(f(x)) ; % determine size of f

h = sqrt(eps) ;
J = zeros(Nf,Nx) ; % initialize the Jacobian

% Compute the Jacobian
for state_index = 1:Nx
        % Create a vector of zeros, with the step size at the
        % appropriate row position
        hvec = zeros(Nx,1) ;
        hvec(state_index) = h ;

        % Compute the gradient of f at x, with the step taken along
        % only one row of x
        f0 = f(x - hvec) ;
        f1 = f(x + hvec) ;
        J(:,state_index) = (f1 - f0)./(2*h) ;
end

end