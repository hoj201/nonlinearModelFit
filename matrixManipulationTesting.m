clear
clc

h = 2 ;
w = 2 ;
d = 3 ;

test = reshape(1:h*w*d,h,w,d)
% reshape(permute(idx,[1 3 2]),[],2)
% test2 = permute(test, [1 3 2]) ;
% reshape(test2,[],w)

test3 = num2cell(test,[1 2])

blkdiag(test3{:})

%% testing equality constraint matrix generation
clear
clc

Ns = 3 ;
Nd = 4 ;
Np = 2 ;
N = Ns*Nd ;

x0 = 5.*ones(Ns,1) ; % fake initial condition
xdata = -2.*ones(Ns,Nd) ; % fake data
x = 2.*ones(Ns,Nd) ; % current system dynamics evaluation

dt = 0.01 ; % fake time delta

dxdt = 3.*ones(Ns,Nd) ;
dfdx = reshape(7.*ones(1,Ns*Ns*Nd),Ns,Ns,Nd) ;
dfdp = reshape(11.*ones(1,Ns*Np*Nd),Ns,Np,Nd) ;

%% EQUALITY CONSTRAINTS
x
dxdt
eq = x - [x0, x(:,1:Nd-1)] - dt.*[zeros(Ns,1),dxdt(:,1:Nd-1)]
% eq = eq(:)'

%% GRADIENT
% reshaping dfdp
dfdpgeq = dt.*reshape(permute(dfdp,[1 3 2]),[],Np) ;

% reshaping dfdx
dfdxcel = num2cell(dfdx,[1 2]) ;
dfdxgeq = dt.*blkdiag(dfdxcel{1:end-1}) ;


geq = [eye(N), [zeros(Ns,Np); dfdpgeq(1:N-Ns,:)]] + ...
      [zeros(Ns,N+Np); -eye(N-Ns)-*dfdxgeq, zeros(N-Ns,Np+Ns)] ;

