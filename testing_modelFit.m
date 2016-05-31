%% Initial testing
clear
clc

f = @testFunction1 ;
p = [2; 4] ;

x0 = [1; -1] ;
p0 = [3; 5] ;
T = 4 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 3 ;

tvec = linspace(0,T,Ntime) ;
uvec = ones(1,Ntime) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
tic
[sol, ~] = user.modelFit() ;
toc
figure(1)
hold on
plot(sol.x')
plot(x')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])
% NOTES: it works, but it's not very accurate


%% Initial testing
clear
clc

f = @testFunction2 ;
p = 10 ;

x0 = -9 ;
p0 = 11 ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 5 ;

tvec = linspace(0,T,Ntime) ;
uvec = ones(1,Ntime) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
[sol, ~] = user.modelFit() ;
plot(x')
hold on
plot(sol.x')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])

% testing dynamics gradient
in = {rand(size(user.x0)),rand(size(user.input(:,1))),rand(size(user.p0))} ;
[~,dfdx,dfdp] = user.fdyn(in{:});
geqAnl = [dfdx,dfdp] ;
geqNum = numericJacobian(@user.fdyn,[1 3],in{:}) ;
norm(geqAnl - geqNum) 

% testing cost gradient
zz = rand(size(user.z0)) ;
[~,geqAnl] = user.cost(zz) ;
geqNum = numericJacobian(@user.cost,1,zz) ;
norm(geqAnl - geqNum)

%% Dependence on number of data points
clear
clc

Ntrials = 5 ;
Ndatalo = 20 ;
Ndatahi = 70 ;

f = @testFunction1 ;
p = [2; 4] ;

x0 = [1; -1] ;
p0 = [3;5] ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
NtimeVec = round(linspace(Ndatalo,Ndatahi,Ntrials)) ;
pvec = zeros(Nparams,Ntrials) ;
ttestvec = zeros(1,Ntrials) ;

figure(2)
hold on

for trial = 1:Ntrials
Ntime = NtimeVec(trial) ;

tvec = linspace(0,T,Ntime) ;
uvec = ones(1,Ntime) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
tic ;
[sol, ~] = user.modelFit() ;
ttestvec(trial) = toc ;
plot(sol.x')
pvec(:,trial) = sol.p ;
end

disp(ttestvec)
disp(pvec)

% NOTES: the time grows a *lot* with the number of data points (duh), and
% rapidly causes fmincon to try to call too many function calls

%% Nonlinear parameter 1
clear
clc

f = @testFunction2 ;
p = 10 ;

x0 = -9 ;
p0 = 15 ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 5 ;

tvec = linspace(0,T,Ntime) ;
uvec = ones(1,Ntime) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
[sol, ~] = user.modelFit() ;

figure(3)
hold on
plot(sol.x')
plot(x')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])


%% Nonlinear parameter 2
clear
clc

f = @testFunction3 ;
p = 0.5 ;

x0 = 5 ;
p0 = 0.7 ;
T = 15 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 100 ;

tvec = linspace(0,T,Ntime) ;
uvec = ones(1,Ntime) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;
% [x2, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p0) ;
% plot(x')
% hold on
% plot(x2')

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
[sol, ~] = user.modelFit() ;
figure(4)
hold on
plot(sol.x')
plot(x')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])

% testing dynamics gradient
in = {rand(size(user.x0)),rand(size(user.input(:,1))),rand(size(user.p0))} ;
[~,dfdx,dfdp] = user.fdyn(in{:});
geqAnl = [dfdx,dfdp] ;
geqNum = numericJacobian(@user.fdyn,[1 3],in{:}) ;
norm(geqAnl - geqNum) 

% testing cost gradient
zz = rand(size(user.z0)) ;
[c,geqAnl] = user.cost(zz) ;
geqNum = numericJacobian(@user.cost,1,zz) ;
norm(geqAnl - geqNum)


%% Sinusoidal input
% clear
% clc

f = @testFunction1 ;
p = [2; 4] ;

x0 = [1; -1] ;
p0 = [1.9; 4.1] ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 200 ;

tvec = linspace(0,T,Ntime) ;
uvec = sin(3*tvec) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
tic
[sol, ~] = user.modelFit() ;
toc
figure(5)
hold on
plot(sol.x')
plot(x')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])


%% Noisy data
clear
clc

f = @testFunction3 ;
p = 0.5 ;

x0 = 5 ;
p0 = 0.7 ;
T = 15 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 150 ;

tvec = linspace(0,T,Ntime) ;
uvec = sin(4*tvec) + 1.01 ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;

% make noisy data
% xnoise = x + 2*rand(size(x)) - 1 ; % uniform noise
xnoise = x + 2*randn(size(x)); % Gaussian noise

% model fitter setup
data = xnoise;
input = uvec ;
Ndata = Ntime ; % currently, this works for fitting to only a single trajectory
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
[sol, ~] = user.modelFit() ;

figure(6)
hold on
plot(x')
plot(xnoise')
plot(sol.x')

legend('Original Data','Noisy Data','Nonlinear Fit')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])

%% Piddling around
clear
clc

f = @testFunction4 ;
p = 1 ;

x0 = [2 -3]' ;
p0 = -1 ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 100 ;

tvec = linspace(0,T,Ntime) ;
uvec = sin(4*tvec) + 1.01 ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;
[x2, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p0) ;

% make noisy data
% xnoise = x + 2*rand(size(x)) - 1 ; % uniform noise
xnoise = x + 2*randn(size(x)); % Gaussian noise

% model fitter setup
data = x;
input = uvec ;
Ndata = Ntime ;
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
[sol, ~] = user.modelFit() ;

figure(7)
hold on
plot(x')
plot(xnoise')
plot(sol.x')

legend('Original Data','Noisy Data','Nonlinear Fit')

disp(['Parameter Fit: ',mat2str(sol.p)])
disp(['Gradient Norm: ',num2str(norm(sol.grad))])
disp(['Hessian Norm: ',num2str(norm(sol.hessian))])

%% Dubins
% clear
% clc

f = @testDubins ;
p = [1 1]' ;

x0 = [0 0 0]' ;
p0 = [1.5 1.5]' ;
T = 5 ;

Nstates = length(x0) ;
Nparams = length(p) ;
Ntime = 100 ;

tvec = linspace(0,T,Ntime) ;
% uvec = ones(1,Ntime);
uvec = sin(4*tvec) ;
[x, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p) ;
[x2, ~, ~, ~] = simulateDynamicsWithInput(f,tvec,uvec,x0,p0) ;

% make noisy data
% xnoise = x + 2*rand(size(x)) - 1 ; % uniform noise
xnoise = x + 0.1*randn(size(x)); % Gaussian noise

% model fitter setup
data = xnoise;
input = uvec ;
Ndata = Ntime ;
user = nonlinearModelFit(f,data,input,x0,p0,T,Nstates,Nparams,Ndata) ;

% testing fmincon call
tic
[sol, ~] = user.modelFit() ;
toc

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

%% testing dynamics gradient
in = {rand(size(user.x0)),rand(size(user.input(:,1))),rand(size(user.p0))} ;
[~,dfdx,dfdp] = user.fdyn(in{:});
geqAnl = [dfdx,dfdp] ;
geqNum = numericJacobian(@user.fdyn,[1 3],in{:}) ;
norm(geqAnl - geqNum) 

%% testing cost gradient
zz = rand(size(user.z0)) ;
[c,geqAnl] = user.cost(zz) ;
geqNum = numericJacobian(@user.cost,1,zz) ;
norm(geqAnl - geqNum)
sol.finalCost

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