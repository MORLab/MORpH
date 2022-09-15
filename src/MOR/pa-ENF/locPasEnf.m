function sysp = locPasEnf(sys,varargin)
% LOCPASENF - Enforce passivity of sys by applying local passivity constraints
%
% Syntax:
%   sysp = locPasEnf(sys)
%   sysp = locPasEnf(sys,opts)
%
% Description:
%       sysp = locPasEnf(sys,opts) perturbs matrix C of sys to obtain a
%       passive model sysp by applying local constraints on the Popov
%       function.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss object, containing LTI system
%       *Optional Input Arguments:*
%       - opts:  structure with execution parameters
%           - .maxIter: 	Maximum number of iterations
%                           [{500} / positive integer]
%           - .step:        Factor controlling how much matrix C is perturbed.
%                           [{0.5} / double]
%           - .verbose:     Control console output.
%                           [{0} / 1]
%           - .makePH:    	Construct port-Hamiltonian system.
%                          	[{false} / true]
%           - .w:           Vector containing sampling points.
%                           [{logspace(-15,5,100)} / vector]
%           - .w0:          Vector containing sampling points passed initially
%                           to samPassive for a finer first sampling.
%                           [{logspace(-15,5,5000)} / vector]
%           -.samPassive.*  Other options that will be passed on to the
%                       	used method (Utility/passivity_check/samPassive);
%                         	Please refer to documentation of the respective
%                         	algorithm.
%           - .ss2phs.*:   	Other options that will be passed on to the
%                          	used method (Utility/ss2phs);
%                        	Please refer to documentation of the respective
%                          	algorithm.
%
% Output Arguments:
%       - sysp:    ssRed object, passive LTI system
%
% See also:
%       sampassive, quadprog, demo_locPasEnf
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

%% Dimensions
numVar = numel(C);

%% Enforcement
% Controllability Gramian
Wc = lyap(A,B*B');
Qc = chol(Wc);
Qcinv = eye(size(Qc))/Qc;

% Calculate passivity violations
[p,v] = samPassive(sys,opts.w0,opts.samPassive);

% Construct new coarse frequency vector
wAdd = v.frequencies(v.violations<0);
w = sort([opts.w,wAdd],'ascend');

% Check passivity
if p
    warning('System is already passive.')
    sysp = sys;
    return
end

% Get data
lambda = v.violations;
f = v.frequencies;
rv = v.eigenvectors;

lambda(lambda>0)=0;

% Initialize Copt
Copt = C;
x0 = zeros(numVar,1);

% Preallocate max violations
max_vio = zeros(opts.maxIter,1);
stopcrit = false;

i=1;
while ~p && i < opts.maxIter
    % Compute matrix Z
    z = zeros(numVar,length(lambda));

    % (Eq.10.29 / 10.40)
    for j=1:length(lambda)
        K = (1i*f(j)*eye(size(A))-A)\B;
        z(:,j) = 2*real(kron(Qcinv.'*K*rv(:,j),conj(rv(:,j))));
    end
    Z=z.';

    % Solve quadratic program
    x = quadprog(2*eye(numVar),zeros(numVar,1),-Z,lambda,[],[],[],[],x0,opts.quadprog);

    % If x=0, solve with old value of C
    if all(x==zeros(size(x)))
        x = quadprog(2*eye(numVar),zeros(numVar,1),-Z,lambda,[],[],[],[],ftv(Copt),opts.quadprog);
    end

    % If still x=0 solve with random initialization
    j=0;
    while all(x==zeros(size(x))) && j<100
        x = quadprog(2*eye(numVar),zeros(numVar,1),-Z,lambda,[],[],[],[],rand(size(x)),opts.quadprog);
    end

    Copt = Copt + opts.step*vtf(x,size(C,1),size(C,2))*Qcinv.';
    sysp = ssRed(A,B,Copt,D,E,'locPasEnf',opts);

    % Calculate passivity violations for perturbed system
    [p,v] = samPassive(sysp,w,opts.samPassive);

    max_vio(i) = min(v.violations);
    if i>5 && abs(max_vio(i)-max_vio(i-5))<1e-10
        stopcrit = true;
        break;
    elseif i>5 && abs(max_vio(i)-max_vio(i-4))<1e-10
        opts.step = 0.9*opts.step;
    end

    % Print current violation to console
    if opts.verbose
        fprintf('Current maximum violation: %.8d\n',min(v.violations))
    end

    if ~p
        % Get data
        lambda = v.violations;
        f = v.frequencies;
        rv = v.eigenvectors;
        lambda(lambda>0)=0;
    end

    i=i+1;
end

if i == opts.maxIter
    fprintf('\t Algorithm converged. Maximum number of iterations reached.\n')
elseif stopcrit
    fprintf('\t Algorithm converged. Stopping criterion triggered.\n')
else
    fprintf('\t Algorithm converged.\n')
end

% Define output
if opts.makePH
    sysp = ss2phs(sysp, opts.ss2phs);
    sysp.parameters = opts;
end

end

function [sys,opts] = parseInputs(sys,varargin)
if nargin < 2
    opts = struct();
else
    opts = varargin{1};
end

OptsAdmissible.samPassive.plot = {false,true};
OptsAdmissible.samPassive.outputVolume = 'all';
OptsAdmissible.w = logspace(-15,5,100);
OptsAdmissible.w0 = logspace(-15,5,5000);
OptsAdmissible.samPassive.tol = 1e-8;
OptsAdmissible.verbose = {0,1};
OptsAdmissible.step = 0.5;
OptsAdmissible.maxIter = 500;
OptsAdmissible.makePH = {false,true};
OptsAdmissible.ss2phs.method = @locPasEnf;

opts = phsMOR_parseOpts(opts, OptsAdmissible);

opts.quadprog = optimoptions('quadprog','Algorithm','active-set','Display','off');

% Check type of sys
if isa(sys,'phs') || isa(sys,'phsRed')
    error('System is already passive.')
elseif ~(isa(sys,'sss') || isa(sys,'ssRed') || isa(sys,'ss'))
    error('Input system is not of type sss, ssRed or ss.')
end

% Check if sys is dae
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('locPasEnf currently only supports ODE systems.')
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