function sysp = hamPasEnf(sys,varargin)
% HAMPASENF - Enforce passivity of sys by applying Hamiltonian passivity
%             constraints.
%
% Syntax:
%   sysp = hamPasEnf(sys)
%   sysp = hamPasEnf(sys,opts)
%
% Description:
%       sysp = hamPasEnf(sys,opts)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss-object, containing LTI system
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .step:            Factor controlling how much matrix C is perturbed.
%                               [{1} / double]
%           - .verbose:         Control console output.
%                               [{0} / 1]
%           - .maxIter:         Maximum number of iterations.
%                               [{1000} / positive integer]
%           - .alpha:           Control parameter for computation of new locations.
%                               [{0.5} / double]
%           - .makePH:          Construct port-Hamiltonian system.
%                               [{false} / true]
%           - .shamPassive.*:   Other options that will be passed on to the
%                               used method (Utility/passivity_check/shamPassive);
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .ss2phs.*:        Other options that will be passed on to the
%                               used method (Utility/ss2phs);
%                               Please refer to documentation of the respective
%                               algorithm.
%
% Output Arguments:
%       - sysp:    sss object, passive LTI system
%
% See also:
%       shampassive, quadprog, demo_hamPasEnf
%
% References:
%       [1] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: 
%           Theory and Applications. Jan. 2016.
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Maximilian Bonauer, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Input parsing
[sys,opts] = parseInputs(sys,varargin{:});

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
E = sys.E;
sysp = sss(A,B,C,D,E);

%% Check definiteness of D+D'
DDT = D + D.';
d = eig(DDT);
if ~(min(d) > 0)
    error('hamPasEnf works only for systems with D+D''>0.')
end

DDTinv = eye(size(DDT))/DDT;

%% Enforcement
% Controllability Gramian
Wc = lyap(A,B*B');
Qc = chol(Wc);
Qcinv = eye(size(Qc))/Qc;

% Calculate passivity violations
[p,v] = shamPassive(sysp,opts.shamPassive);

% Check passivity
if p
    warning('System is already passive. Returning original system.')
    sysp.D=sys.D;
    return
end

% Get purely imaginary eigenvalues and corresponding eigenvectors
nue = v.eigenvectors;
omega = v.w;
slopes = v.slopes;

% Dimensions
n = size(A,1);
m = size(B,2);

% Initialize Copt
Copt = C;
numVar = numel(C);
x0 = zeros(numel(C),1);

% Preallocate max violations
max_vio = zeros(opts.maxIter,1);
stopcrit = false;

i=1;
while ~p && i < opts.maxIter
    % Preallocation
    z = zeros(n*m,length(omega));
    eta = zeros(1,length(omega));

    % Construct optimization problem
    for j=1:length(omega)
        Xr = nue(:,j);

        % Compute target eigenvalue
        if slopes(j) > 0
            if j==1
                omega_new = omega(j) - opts.alpha*(omega(j)-0);
            else
                omega_new = omega(j) - opts.alpha*(omega(j)-omega(j-1));
            end
        else
            if j==length(omega)
                omega_new = omega(j) + opts.alpha*(2*omega(j)-omega(j));
            else
                omega_new = omega(j) + opts.alpha*(omega(j+1)-omega(j));
            end
        end

        % Compute constraints
        y = DDTinv*(Copt*Xr(1:n) + B.'*Xr(n+1:end));
        z(:,j) = real(kron(Qcinv.'*Xr(1:n),conj(y)));
        eta(:,j) = (omega(j)-omega_new)*imag(Xr(1:n)'*Xr(n+1:end));
    end

    % Compute optimal solution
    if cond(z.'*z) < 10 % Tuneable Parameter
        x = z/(z.'*z)*eta.';
    else
        x = quadprog(2*eye(numVar),zeros(numVar,1),[],[],z',eta,[],[],x0,opts.quadprog);
        if all(x==zeros(size(x)))
            x = quadprog(2*eye(numVar),zeros(numVar,1),[],[],z',eta,[],[],ftv(Copt),opts.quadprog);
        end
        j=0;
        while all(x==zeros(size(x))) && j<100
            x = quadprog(2*eye(numVar),zeros(numVar,1),[],[],z',eta,[],[],rand(size(x)),opts.quadprog);
        end
    end

    Copt = Copt + opts.step*vtf(x,size(C,1),size(C,2))*Qcinv.';
    sysp = ssRed(A,B,Copt,D,E,'hamPasEnf',opts);

    % Calculate passivity violations for perturbed system
    [p,v] = shamPassive(sysp,opts.shamPassive);

    % Check convergence
    if ~p
        max_vio(i) = min(v.violations);
        if i>5 && (abs(max_vio(i)-max_vio(i-5))<1e-10 || abs(max_vio(i)-max_vio(i-5))>1e10)
            stopcrit = true;
            break;
        end
    end

    % Disp current violations
    if opts.verbose
        fprintf('Current maximum violation: %d.8\n',min(v.violations))
    end

    if ~p
        % Update data
        nue = v.eigenvectors;
        omega = v.w;
        slopes = v.slopes;
    end

    i=i+1;
end

% Print convergence info to console
if i == opts.maxIter
    fprintf('\t Algorithm converged. Maximum number of iterations reached.\n')
elseif stopcrit
    fprintf('\t Algorithm converged. Stopping criterion triggered.\n')
else
    fprintf('\t Algorithm converged.\n')
end

% Remove artificial feedthrough
if ~(min(d) > 0)
    sysp.D = zeros(size(D));
end

% Define output
if opts.makePH
    sysp = ss2phs(sysp, opts.ss2phs);
    %     sysr.method = @hamPasEnf;
    sysp.parameters = opts;
end

end

%% Auxilliary function
function [sys,opts] = parseInputs(sys,varargin)
if nargin < 2
    opts = struct();
else
    opts = varargin{1};
end

% Default values
OptsAdmissible.shamPassive = struct();
OptsAdmissible.step = 1;
OptsAdmissible.verbose = {0,1};
OptsAdmissible.maxIter = 1000;
OptsAdmissible.alpha = 0.5;
OptsAdmissible.makePH = {false,true};
OptsAdmissible.ss2phs = struct();

opts = phsMOR_parseOpts(opts, OptsAdmissible);

% Quadprog options
opts.quadprog = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

% Check type of sys
if isa(sys,'phs') || isa(sys,'phsRed')
    error('System is already passive.')
elseif ~(isa(sys,'sss') || isa(sys,'ssRed') || isa(sys,'ss'))
    error('Input system is not of type sss, ssRed or ss.')
end

% Check if sys is dae
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('hamPasEnf currently only supports ODE systems.')
end

% Check if E exists
if ~isempty(sys.E)
    % Make system explicit
    sys.A = sys.E\sys.A;
    sys.B = sys.E\sys.B;
    sys.E = eye(size(sys.A));
else
    sys.E = eye(size(sys.A));
end

end
