function sysp = prlPasEnf(sys,varargin)
% PRLPASENF - Enforce passivity of LTI system by passivity enforcement via
%             Positive Real Lemma.
%
% Syntax:
%   sysp = prlPasEnf(sys)
%   sysp = prlPasEnf(sys, opts)
%
% Description:
%       sysp = prlPasEnf(sys, opts) perturbs matrix C of sys to obtain a
%       passive model sysp by searching for the solution of the modified
%       KYP inequality.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss object, containing LTI system
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .tol:            Tolerance for solution X to KYP inequality to be
%                               positive definite.
%                               [{1e-8} / positive double]
%           - .verbosity:       Regulates ouput of cvx in command window.
%                               [{0} / 1]
%           - .solver:          Define the solver for the Positive Real lemma.
%                               [{'cvx'} / 'yalmip']
%           - .makePH:          Construct port-Hamiltonian system.
%                               [{false} / true]
%           - .ss2phs.*:        Other options that will be passed on to the
%                               used method (Utility/ss2phs);
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .yalmip.solver:   Solver option that will be passed on to sdpsettings
%                               to solve the PRL LMI.
%                               Please refer to documentation of the respective
%                               algorithm.
%                               [{'sedumi'} / string]
%
%
% Output Arguments:
%       - sysp:    ssRed object, passive LTI system
%
% See also:
%       cvx, yalmip, demo_prlPasEnf
%
% References:
%       [1] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: 
%           Theory and Applications. Jan. 2016.
%       [2] M. Grant and S. Boyd. CVX: Matlab Software for Disciplined Convex 
%           Programming, version 2.1. Mar. 2014. url: http://cvxr.com/cvx.
%       [3] M. Grant and S. Boyd. “Graph implementations for
%           nonsmooth convex programs." In: Recent Advances in
%           Learning and Control. Ed. by V. Blondel, S. Boyd, and H.
%           Kimura. Lecture Notes in Control and Information Sciences.
%           Springer-Verlag Limited, 2008, pp. 95–110.
%       [4] J. Löfberg. “YALMIP : A Toolbox for Modeling and Optimization in MATLAB.” 
%           In: In Proceedings of the CACSD Conference. Taipei, Taiwan, 2004.
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
%               If you use the software packages 'CVX' [2,3] or 'YALMIP' [4], 
%               please refer to their specific license.
%-----------------------------------------------------------------------

%% Input parsing
[sys,opts] = parseInputs(sys,varargin{:});

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

%% Dimensions
n = size(A,1);
m = size(B,2);

%% Passivity Enforcement
% Controllability Gramian
Wc = lyap(A,B*B');
try
    Qc = chol(Wc);
catch
    error('System is not controllable and/or stable.')
end

Qcinv = eye(size(Qc))/Qc;
% Solve semidefinite program
if strcmp(opts.solver,'cvx')
    if opts.verbosity==1
        cvx_begin sdp
    else
        cvx_begin quiet sdp
    end
    %         cvx_precision high
    variable X(n,n);
    variable Xi(m,n);
    minimize(sum_square(vec(Xi)));
    subject to
    -[A'*X+X'*A, X*B-C'-Qcinv*Xi'; B'*X-C-Xi*Qcinv', -D-D'] == semidefinite(n+m);
    X - opts.tol*eye(n) == semidefinite(n);
    cvx_end

    % Give feedback on solution accuracy
    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        error('CVX failed to compute a solution to the PRL.')
    elseif strcmp(cvx_status,'Inaccurate/Solved')
        warning('CVX solution to the PRL is inaccurate.')
    end
else
    yalmip('clear')
    opt = sdpsettings('solver',opts.yalmip.solver,'verbose',opts.verbosity); % Optimization settings

    % Define variables
    X = sdpvar(n,n);
    Xi = sdpvar(m,n);

    % Define matrices
    mat = [A'*X+X'*A, X*B-C'-Qcinv*Xi'; B'*X-C-Xi*Qcinv', -D-D'];

    % Define Constraints
    F = [];
    F = [F X>=opts.tol*eye(n)];
    F = [F mat<=-opts.tol*eye(size(mat))];

    % Optimization Problem
    info = optimize(F,norm(Xi,'fro')^2,opt); 
    
    Xi = value(Xi);
    X = value(X);
    
    if info.problem == 0
        % Nothing to do
    elseif info.problem == 1
        warning('YALMIP thinks problem is infeasible.')
    else
        warning('Solution to KYP inequality may be inaccurate.')
    end
    
    if any(isnan(X))
        error('No accurate solution to KYP inequality found.')
    end
end

% Add perturbation to C
C = C+Xi*Qcinv';

% Output
opts.X = X;

% Create ss object
sysp = ssRed(A,B,C,D,eye(size(A)),'prlPasEnf',opts);

% Define output
if opts.makePH
    sysp = ss2phs(sysp, X, opts.ss2phs);
    sysp.parameters = opts;
end

end

function [sys,opts] = parseInputs(sys,varargin)
if nargin < 2
    opts = struct();
else
    opts = varargin{1};
end

% Default options
OptsAdmissible.verbosity = {0,1};
OptsAdmissible.tol = 1e-8;
OptsAdmissible.makePH = {false,true};
OptsAdmissible.ss2phs.method = @prlPasEnf;
OptsAdmissible.solver = {'cvx','yalmip'};
OptsAdmissible.yalmip.solver = 'sedumi';

opts = phsMOR_parseOpts(opts, OptsAdmissible);

% Check type of sys
if isa(sys,'phs') || isa(sys,'phsRed')
    error('System is already passive.')
elseif ~(isa(sys,'sss') || isa(sys,'ssRed') || isa(sys,'ss'))
    error('Input system is not of type sss, ssRed or ss.')
end

% Check if sys is DAE
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('prlPasEnf currently only supports ODE systems.')
end

% Check if E exists
if ~isempty(sys.E)
    % Make system explicit
    sys.A = sys.E\sys.A;
    sys.B = sys.E\sys.B;
    sys.E = eye(size(sys.A));
end

% Check third party software
switch opts.solver
    case 'cvx'
        thirdPartyCheck('CVX');
    case 'yalmip'
        thirdPartyCheck('YALMIP');
end

end