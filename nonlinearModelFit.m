classdef nonlinearModelFit
% TO DO
% - maybe the data instances should always be in the 3rd dimension, instead
%   of the 2nd dimension for *just* x and dxdt
% - add ability to change options
% - enable varying time steps?

   properties
        % The user enters the following properties into the config file:
        fdyn      % system dynamics and derivatives wrt states and params (function)
        Nstates   % number of states of system (int)
        Nparams   % number of unknown parameters to identify (int)
        Ndata     % number of data points (integer)
        data      % state data matrix (Nstates x Ndata)
        input     % input data matrix (Ninput x Ndata)
        x0        % initial condition (should correspond to first row of data)
        p0        % initial guess at parameters
        T         % time span of data

        % These properties are created and used internally
        z0        % initial condition of decision variable z
        dt        % time step between data points (T/Ndata)
   end
   methods
    % CONSTRUCTOR
        function user = nonlinearModelFit(dynamics_function,data,input,...
                            x0,p0,T,Nstates,Nparams,Ndata)
           % Save user input
           user.fdyn = dynamics_function ;
           user.Nstates = Nstates ;
           user.Nparams = Nparams ;
           user.data = data ;
           user.input = input ;
           user.Ndata = Ndata ;
           user.x0 = x0 ;
           user.p0 = p0 ;
           user.T = T ;
           
           % Create useful variables
           user.dt = T/Ndata ;
           [x, ~, ~, ~] = user.forwardSimulateDynamics(user.p0) ;
           user.z0 = [x(:); p0] ;
        end
        
    % OPTIMIZATION CALL
        function [solution, problem] = modelFit(user)
            options = optimoptions('fmincon','GradObj','on') ;
            options.MaxFunEvals = 100000 ;
            problem.objective = @user.cost ;
            problem.x0 = user.z0 ;
            problem.nonlcon = @user.nonlinearConstraints ;
            problem.solver = 'fmincon' ;
            problem.options = options ;
            
            [z,fval,exitflag,output,lambda,grad,hessian] = fmincon(problem) ;
            [x,p] = user.statedecode(z) ;
            solution.x = x ;
            solution.p = p ;
            solution.finalCost = fval ;
            solution.exitflag = exitflag ;
            solution.output = output ;
            solution.lambda = lambda ;
            solution.grad = grad ;
            solution.hessian = hessian ;
        end

    % METHODS TO USE IN OPTIMIZATION
        function [c,gc] = cost(user,z)
            % Calculate the cost associated with the decision variable vector
            % z, along with the gradient of the cost, at a single optimization
            % step. The vector z is a column vector of length (Ns*Nd + Np)

            % Extract states (Nstates x Ndata) and parameters (Nparams x 1)
            [x,~] = user.statedecode(z) ;
            c = (sum((x(:)-user.data(:)).^2)) ;
            gc = [2.*(x(:)-user.data(:))', zeros(1,user.Nparams)] ;
        end

        function [cons,eq,gcons,geq] = nonlinearConstraints(user,z)
        % Given the decision variable vector z, set up the optimization
        % constraints and the gradient of these constraints

        % Set inequality constraints to null
            cons = [] ;
            gcons = [] ;
        % Get equality constraints (this is a separate function, so that
        % the user can test the accuracy of the analytic gradients that
        % should be provided in user.fdyn)
            [eq, geq] = user.equalityConstraints(z) ;
            eq = eq' ;
            geq = geq' ;
        end
        
        function [eq,geq] = equalityConstraints(user,z)
        % Given the decision variable vector z, calculate the equality
        % constraints and the equality constraint gradient using the
        % user-provided system dynamics
        
        % Extract states (Nstates x Ndata) and parameters (Nparams x 1)
            Ns = user.Nstates ;
            Np = user.Nparams ;
            Nd = user.Ndata ;
            N = Ns*Nd ;
            [x,p] = user.statedecode(z) ;
            
        % Evaluate the system dynamics and get the relevant gradients
            [dxdt, dfdx, dfdp] = user.evaluateDynamics(x,p);

        % Calculate the equality constraints
            eq = x - [user.x0, x(:,1:Nd-1)] - ...
                 user.dt.*[zeros(Ns,1),dxdt(:,1:Nd-1)] ;
            
            % single point with dynamics
%             eq = x - [user.x0, x(:,1:Nd-1)] ;
%             datumIdx = 2 ;
%             rowIdx = datumIdx*Ns+1 ;

            % dynamics equal to zero
%             eq(rowIdx) = eq(rowIdx) - user.dt.*dxdt(1,datumIdx) ;
%             eq = -user.dt.*[zeros(Ns,1),dxdt(:,1:Nd-1)] ;
            eq = eq(:) ;
            
        % Calculate the equality constraint gradient
            % Reshape dfdx
            dfdxcel = num2cell(dfdx,[1 2]) ;
            dfdxgeq = -user.dt.*blkdiag(dfdxcel{1:end-1}) ;
            dfdxgeq2 = -user.dt.*blkdiag(dfdxcel{2:end-1}) ;
            
            % Reshape dfdp
            dfdpgeq = -user.dt.*reshape(permute(dfdp,[1 3 2]),[],Np) ;

            % Generate equality constraint gradient
            % boring little one
            geq = [eye(N), [zeros(Ns,Np); dfdpgeq(1:N-Ns,:)]] + ...
                  [zeros(Ns,N+Np); -eye(N-Ns)+dfdxgeq, zeros(N-Ns,Np+Ns)] ;

            % pick just a single eq constraint!
%             geq = [eye(N), zeros(N,Np)] + ...
%                   [zeros(Ns,N+Np); -eye(N-Ns), zeros(N-Ns,Np+Ns)] ;% + ...
%             geq(rowIdx) = geq(rowIdx) - user.dt.*dfdx(1,1,datumIdx) ;
%             geq(rowIdx,N+1) = geq(rowIdx,N+1) -user.dt.*dfdp(1,1,datumIdx) ;

            % the big ugly one that works
%             geq = [eye(N), [zeros(Ns,Np);dfdpgeq(1:N-Ns,:)]] + ...
%                   [zeros(Ns,N+Np); -eye(N-Ns)+dfdxgeq, zeros(N-Ns,Np+Ns)] + ...
%                   [zeros(2*Ns,N+Np); dfdxgeq2, zeros(N-2*Ns,Np+2*Ns)] + ...
%                   [zeros(2*Ns,N+Np); zeros(N-2*Ns,Ns), -dfdxgeq2, zeros(N-2*Ns,Np+Ns)] ;
        end

    % UTILITY FUNCTIONS
        function [dxdt, dfdx, dfdp] = evaluateDynamics(user,x,p)
        % Given the system state x and parameters p, evaluate the system
        % dynamics and return dx/dt and the gradients of the dynamics
        % with respect to x and p (df/dx and df/dp)
            u = user.input ;
            Ns = user.Nstates ;
            Np = user.Nparams ;
            Nd = user.Ndata ;
            
            % set up matrices (will need to add an extra dimension when the
            % time comes to set up multiple input observations with Ntime
            % and Ndata as the two data dimensions)
            dxdt = zeros(Ns,Nd) ;
            dfdx = zeros(Ns,Ns,Nd) ;
            dfdp = zeros(Ns,Np,Nd) ;
            
            % evaluate at initial conditions
            [dxdtPoint,dfdxPoint,dfdpPoint] = user.fdyn(x(1:Ns),u(:,1),p) ;
            dxdt(:,1) = dxdtPoint ;
            dfdx(:,:,1) = dfdxPoint ;
            dfdp(:,:,1) = dfdpPoint ;
            
            % this can be improved with arrayfun probably
            for idx = 2:Nd
                % evaluate the system dynamics at each data point
                [dxdtPoint,dfdxPoint,dfdpPoint] = user.fdyn(x(:,idx-1),u(:,idx),p) ;
                
                % fill in the matrices!
                dxdt(:,idx) = dxdtPoint ;
                dfdx(:,:,idx) = dfdxPoint ;
                dfdp(:,:,idx) = dfdpPoint ;
            end
        end
        
        function [x, dxdt, dfdx, dfdp] = forwardSimulateDynamics(user,p)
        % THIS IS PRETTY MUCH THE SAME AS EVALUATEDYNAMICS
        % Given some system parameters p, simulate the system dynamics and
        % return the trajectory, dx/dt, and the gradients of the dynamics
        % with respect to x and p (df/dx and df/dp)
            u = user.input ;
            Ns = user.Nstates ;
            Np = user.Nparams ;
            Nd = user.Ndata ;
            
            % set up matrices (will need to add an extra dimension when the
            % time comes to set up multiple input observations with Ntime
            % and Ndata as the two data dimensions)
            x = [user.x0, zeros(Ns,Nd-1)] ;
            dxdt = zeros(Ns,Nd) ;
            dfdx = zeros(Ns,Ns,Nd) ;
            dfdp = zeros(Ns,Np,Nd) ;
            
            % evaluate at initial conditions
            [dxdtPoint,dfdxPoint,dfdpPoint] = user.fdyn(x(1:Ns),u(:,1),p) ;
            dxdt(:,1) = dxdtPoint ;
            dfdx(:,:,1) = dfdxPoint ;
            dfdp(:,:,1) = dfdpPoint ;
            
            % this can be improved with arrayfun probably
            for idx = 2:Nd
                % evaluate the system dynamics at each data point
                [dxdtPoint,dfdxPoint,dfdpPoint] = user.fdyn(x(:,idx-1),u(:,idx-1),p) ;
                
                % fill in the matrices!
                x(:,idx) = x(:,idx-1) + user.dt.*dxdtPoint ;
                dxdt(:,idx) = dxdtPoint ;
                dfdx(:,:,idx) = dfdxPoint ;
                dfdp(:,:,idx) = dfdpPoint ;
            end
        end

        function [x,p] = statedecode(user,z)
        % Given the decision variables vector z and a datum index n (in 1 to
        % Ndata), return the state of the system and the current parameters
           Ns = user.Nstates ;
           Np = user.Nparams ;
           Nd = user.Ndata ;
           x = reshape(z(1:Ns*Nd),Ns,Nd) ;
           p = z(end-Np+1:end) ;
        end

        function z = stateencode(~,x,p)
        % Given the state at and the current parameters p, return an
        % updated decision variable vector z
           z = [x(:);p] ;
        end
    end
end