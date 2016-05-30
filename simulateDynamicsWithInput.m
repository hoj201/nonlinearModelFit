function [x, dxdt, dfdx, dfdp] = simulateDynamicsWithInput(fxup,tvec,uvec,x0,p)
% Given a vector function f(x,u,p), a vector of time instances, a vector
% of inputs corresponding to those time instances, an initial condition,
% and system parameters, forward-simulate the function from some initial
% condition
%
% INPUTS:
%   fxup    -   system dynamics, which takes in the state x, input u, and
%               system parameters p, and returns dx/dt
%   tvec    -   time vector of artbitrary length Ntime, with the first
%               element assumed to be 0
%   uvec    -   vector of system inputs at each time instance of tvec,
%               assumed to be of size (Ninputs x Ttime)
%   x0      -   initial condition, of size (Nstates x 1)
%   p       -   system parameters, of size (Nparams x 1)
%
% OUTPUTS:
%   x       -   system trajectory, of size (Nstates x Ntime)
%   xd      -   dx/dt, of size (Nstates x Ntime)

    T = length(tvec) ; % number of time steps to simulate
    Ns = length(x0) ; % number of states of system
    Np = length(p) ; % number of unknown parameters
    x = [x0, zeros(Ns,T-1)] ;
    dxdt = zeros(Ns,T) ;
    dfdx = zeros(Ns,Ns,T) ;
    dfdp = zeros(Ns,Np,T) ;
    [dxdt0, dfdx0, dfdp0] = fxup(x0,uvec(:,1),p) ;
    dxdt(:,1) = dxdt0 ;
    dfdx(:,:,1) = dfdx0 ;
    dfdp(:,:,1) = dfdp0 ;
    for dt = 2:T
        xold = x(:,dt-1) ;
        tstep = tvec(dt) - tvec(dt-1) ;
        [dxdtPoint, dfdxPoint, dfdpPoint] = fxup(xold,uvec(:,dt-1),p) ;
        dxdt(:,dt) = dxdtPoint ;
        dfdx(:,:,dt) = dfdxPoint ;
        dfdp(:,:,dt) = dfdpPoint ;
        x(:,dt) = xold + tstep*dxdtPoint ;
    end
end