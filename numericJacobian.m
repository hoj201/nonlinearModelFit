function J = numericJacobian(f,dfIdx,varargin)
% Evaluate the numeric Jacobian of a function f, which can take Nfin inputs
% of varying sizes, and returns a vector of size Nfout x 1.
%
% INPUTS:
%   f        Nfout x 1 vector function of form f(x1,x2,...), where each xi
%            input is a vector sized appropriately for f.
%
%   dfIdx    A vector of indices expressing which inputs of f to take
%            the Jacobian with respect to. For example, if dfIdx is [1 3],
%            then the Jacobian of f will be evaluated for only the 1st and
%            3rd arguments to f, and will be of dimension Nfout x (Nx1 +
%            Nx3), where Nx1 and Nx3 are the lengths of the 1st and 3rd
%            argument, respectively.
%
%   x1,x2... The inputs to f, expressing the point x = {x1,x2,...} about
%            which to evaluate jacobian(f(x1,x2,...)). These should be
%            passed to numericJacobian as they would be passed to f.
%
% OUTPUTS:
%   J        The Jacobian of f evaluated at the point x = {x1,x2,...}, with
%            respect to the inputs listed by dfIdx.

% Determine sizes of inputs
numInputsToEval = length(dfIdx) ;
[dfInputs{1:numInputsToEval}] = varargin{dfIdx} ;
Nxi = cellfun('length',dfInputs) ; % length of each input ro f
Nxtotal = sum(Nxi) ; % total number of columns of J
Nf = length(f(varargin{:})) ; % total number of rows of J
columnOffsets = [0 Nxi(1:end-1)] ;

% Initialize step size and empty Jacobian
h = sqrt(eps) ;
J = zeros(Nf,Nxtotal) ;

for inputCount = 1:numInputsToEval
    Nx = Nxi(inputCount) ;
    inputIndex = dfIdx(inputCount) ;
    colOff = sum(columnOffsets(1:inputCount)) ;
    hmat = diag(h.*ones(Nx,1)) ;
    for stateIndex = 1:Nx
        % Create a vector of zeros, with the step size at the
        % appropriate row position
        hvec = hmat(:,stateIndex) ;
        
        % Construct inputs to f with step taken only along the particular
        % parameter corresponding to the current column of J
        fInMinus = varargin ;
        fInMinus{inputIndex} = varargin{inputIndex} - hvec ;
        
        fInPlus = varargin ;
        fInPlus{inputIndex} = varargin{inputIndex} + hvec ;
        
        % Compute the gradient of f at xi, with the step taken along
        % only one row of xi
        fMinus = f(fInMinus{:}) ;
        fPlus = f(fInPlus{:}) ;
        
        % Fill in J
        J(:,stateIndex + colOff) = (fPlus - fMinus)./(2*h) ;
    end
end
        