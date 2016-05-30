f = @dubinsCarDynamics ;
p = [1 1 1]' ;

x0 = [0 0 0]' ;
p0 = [1.5 1.5 1.5]' ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 75 ;

tvec = linspace(0,T,Ntime) ;
uvec = [5*cos(10.*tvec); ones(1,Ntime)];
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;
[x2, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p0) ;

%% add noise to data
xnoise = x + 0.25*(rand(size(x)) - 0.5) ; % uniform noise
% xnoise = x + 0.01*randn(size(x)); % Gaussian noise

% model fitter setup
data = xnoise;
input = uvec ;
Ndata = Ntime ;
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

%% testing dynamics gradient
in = {rand(size(user.x0)),rand(size(user.input(:,1))),rand(size(user.p0))} ;
[~,dfdx,dfdp] = user.fdyn(in{:});
geqAnl = [dfdx,dfdp] ;
geqNum = numericJacobian(@user.fdyn,[1 3],in{:}) ;
norm(geqAnl - geqNum) 

%% testing cost gradient
zz = rand(size(user.z0)) ;
[~,geqAnl] = user.cost(zz) ;
geqNum = numericJacobian(@user.cost,1,zz) ;
norm(geqAnl - geqNum)
% norm(geqAnl - geqNum)*sol.finalCost

%% testing equality gradient
u = user.input ;
Ns = user.Nstates ;
Np = user.Nparams ;
Nd = user.Ndata ;

[r,c] = size(user.z0) ;
v = randn(r,c) ;
[eq,geqAnl] = user.equalityConstraints(v) ;
tic
geqNum = numericJacobian(@user.equalityConstraints,1,v) ;
toc
disp(['Numeric vs. Analytic Gradients: ', num2str(norm(geqAnl - geqNum))])

%% parameter fitting
tic
[sol, ~] = user.modelFit() ;
toc

%%
figure(8)
hold on
plot(x(1,:),x(2,:))
plot(xnoise(1,:),xnoise(2,:),':')
plot(sol.x(1,:),sol.x(2,:))
axis square

legend('Original Data','Noisy Data','Fit Model')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])