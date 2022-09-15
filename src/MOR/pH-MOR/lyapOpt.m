function [redSys] = lyapOpt(sys,redOrder, varargin)
% LYAPOPT - Optimizes the matrices J, R and G of a
% reduced-order model on the product manifold M = sym+(r) x skew(r) x
% R^{rxm} in an H2-optimal way using the Lyapunov framework.
%
% Syntax:
%   redSys = LYAPOPT(sys,redOrder)
%   redSys = LYAPOPT(sys,redOrder,redSys0)
%   redSys = LYAPOPT(sys,redOrder,redSys0,Opts)
%
% Description:
%       redSys = lyapOpt(sys, redOrder) returns a reduced PH system of order
%       redOrder.
%
%       Argument Opts can be used for further adjustments.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object (large-scale model)
%       - redOrder:	desired order of the reduced-order model
%
%       *Optional Input Arguments:*
%       - redSys0: Initial phs system used as a starting point on the manifold
%       - Opts:  structure with execution parameters
%           - checkGradient:    Numerical check of gradient
%                               [{false} / true]
%           - checkHessian:     Numerical check of Hessian
%                               [{false} / true]
%           - lyapSolver:       Solver used for Lyapunov equations
%                               [{'mess'} / 'lyap']
%           - .manopt.*         Other options that will be passed on to the 
%                               third-party software 'Manopt'
%                               Please refer to its documentation (doc manopt).
%           - .phs.*            Other options that will be passed on to the class 'phs';
%                               Please refer to its documentation (doc phs).
%
% Output Arguments:
%       - redSys: 	phsRed object (reduced-order model)
%
% See Also:
%       prOpt
%
% References:
%       [1] K. Sato. “Riemannian Optimal Model Reduction of Linear
%           Port-Hamiltonian Systems.” In: Automatica J. IFAC 93
%           (2018), pp. 428–434.
%       [2] N. Boumal et al. “Manopt, a Matlab Toolbox for Optimization on Manifolds.” 
%           In: Journal of Machine Learning Research 15.42 (2014), pp. 1455–1459. 
%           url: https://www.manopt.org.
%       [3] J. Saak, M. Köhler, and P. Benner. M-M.E.S.S.-2.1 – The
%           Matrix Equations Sparse Solvers library.
%           URL: https://www.mpi-magdeburg.mpg.de/projects/mess
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%               This function uses the software packages 'Manopt' [2] and 'M-M.E.S.S.' [3]; 
%               please refer to their specific license.
%-----------------------------------------------------------------------

%% Input parsing
[redSys0, Opts] = parseInputs(sys, redOrder, varargin{:});

%% Transform original system to explicit form (E=I)
if sys.isImplicit
    warning('off', 'phs:phs:ChangedProperty');
    sys = makeExplicit(sys);
    warning('on', 'phs:phs:ChangedProperty');
end

%% Initialization
% Create the problem structure
problem.M = productmanifold(struct('Jr', skewsymmetricfactory(redOrder),...
    'Rr', sympositivedefinitefactory(redOrder),...
    'Gr', euclideanfactory(redOrder,size(sys.G,2))));

if isa(redSys0,'phs')
    if redSys0.isImplicit || ~isequal(redSys0.Q,speye(redSys0.dim))
        % Transform redSys0 to E=Q=I
        warning('off', 'phs:phs:ChangedProperty');
        redSys0 = scaling(redSys0);
        redSys0 = makeExplicit(redSys0,struct('keepScaled',true));
        warning('on', 'phs:phs:ChangedProperty');
    end
    [x0.Jr, x0.Rr, ~, x0.Gr] = getMatrices(redSys0);

    % Enforce positive definiteness of x0.Rr
    x0.Rr = makePSD(x0.Rr,1e-12);

else
    % random initialization
    x0 = problem.M.rand();
end

%% Link cost function, gradient and Hessian
problem.cost = @(x)cost(x,sys,redOrder,Opts);
problem.grad = @(x)grad(x,sys,redOrder,Opts);
problem.M.inner = @metric;
problem.hess = @(x,u)hess(x,u,sys,redOrder,Opts);

%% Check consistency of gradient and Hessian
if Opts.checkGradient
    checkgradient(problem);
end
if Opts.checkHessian
    checkhessian(problem);
end

%% Trust region algorithm
if Opts.manopt.maxiter > 0
    [x, ~, info] = trustregions(problem,x0,Opts.manopt);
else
    warning('MORpH:lyapOpt:wrongInput',['The desired maximum number of iterations is smaller than 1.\n',...
        'Exit with the initial model.']);
    x = x0;
    info = [];
end

%% Create redSys object
redSys = phsRed(x.Jr, x.Rr, eye(redOrder), x.Gr,Opts.phs);
redSys.method = @lyapOpt;
redSys.parameters = Opts;
redSys.info = info;

end

function [redSys0,Opts] = parseInputs(sys, redOrder, varargin)

% Check number of inputs provided
narginchk(2,4);

% Check original system
if ~isa(sys,'phs') || sys.isDAE || ...
        ~isequal(sys.P,zeros(size(sys.P))) || ~isequal(sys.S,zeros(size(sys.S))) || ~isequal(sys.N,zeros(size(sys.N)))
    error('MORpH:lyapOpt:wrongInput', ...
        ['This algorithm currently only works for pHODEs of the following form:\n',...
        'E*dx/dt = (J - R)*Q*x(t) + G*u(t)\n',...
        'y(t) = G^T*Q*x(t)']);
end

% Check redOrder
if redOrder > sys.dim
    error('MORpH:lyapOpt:wrongInput', 'Reduced Order should be smaller than original system order')
end

% Option parsing
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.phs = struct;
OptsAdmissible.checkGradient = {false,true};
OptsAdmissible.checkHessian = {false,true};
OptsAdmissible.lyapSolver = {'mess','lyap'};
OptsAdmissible.manopt.tolgradnorm = 1e-8;
OptsAdmissible.manopt.maxiter = 1000;
% OptsAdmissible.manopt.maxinner = 50;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% Parse optional arguments
if length(varargin)==1
    redSys0 = varargin{1};
else
    redSys0 = false;
end

% Some checks if an initial reduced model is provided
if isobject(redSys0)
    if ~isa(redSys0,'phs') || redSys0.isDAE || ...
            ~isequal(redSys0.P,zeros(size(redSys0.P))) || ~isequal(redSys0.S,zeros(size(redSys0.S))) || ~isequal(redSys0.N,zeros(size(redSys0.N)))
        error('MORpH:lyapOpt:wrongInput', ...
            ['Initial model must be a pHODE of the following form:\n',...
            'E*dx/dt = (J - R)*Q*x(t) + G*u(t)\n',...
            'y = G^T*Q*x']);
    elseif redSys0.dim ~= redOrder || size(redSys0.G,2)~=size(sys.G,2)
        error('MORpH:lyapOpt:wrongInput', ...
            ['Initial model should have the dimension redOrder and the same number of inputs as the original model.']);
    end
end

% Check for 3rd-party software
thirdPartyCheck('Manopt');
if strcmp(Opts.lyapSolver,'mess')
    thirdPartyCheck('M.E.S.S.');
end

end

function [metr] = metric(x,u,v)

metr = trace(u.Jr'*v.Jr) ...
    + trace((x.Rr\u.Rr)*(x.Rr\v.Rr)) ...
    + trace(u.Gr'*v.Gr);

end

function [f] = cost(x,sys,redOrder,Opts)

Jr = x.Jr; Rr = x.Rr; Gr = x.Gr;
A = (sys.J-sys.R)*sys.Q; G = sys.G;

if strcmp(Opts.lyapSolver,'mess')
    % Error system
    Ae = [A zeros(sys.dim,redOrder);zeros(redOrder,sys.dim) Jr-Rr];
    Be = [G;Gr];
    % Solve Lyapunov equations of error system
    Z = mess_lyap(Ae, Be); Ec = Z*Z';
    X = Ec(1:sys.dim,sys.dim+1:end);
    P = Ec(sys.dim+1:end,sys.dim+1:end);
else
    P = lyap(Jr-Rr,(Jr-Rr)',Gr*Gr');
    X = lyap(A,(Jr-Rr)',G*Gr');
end

f = trace(Gr'*P*Gr - 2*Gr'*X'*G);

end

function [grad] = grad(x,sys,redOrder,Opts)

Jr = x.Jr; Rr = x.Rr; Gr = x.Gr;
A = (sys.J-sys.R)*sys.Q; G = sys.G;

% Solve Lyapunov equations
if strcmp(Opts.lyapSolver,'mess') % MESS solver

    % Error system
    Ae = [A zeros(sys.dim,redOrder);zeros(redOrder,sys.dim) Jr-Rr];
    Be = [G;Gr];
    Ce = [G' -Gr'];

    % Solve Lyapunov equations of error system
    Z = mess_lyap(Ae, Be); Ec = Z*Z';
    Z = mess_lyap(Ae', Ce'); Eo = Z*Z';
    P = Ec(sys.dim+1:end,sys.dim+1:end);
    Q = Eo(sys.dim+1:end,sys.dim+1:end);
    X = Ec(1:sys.dim,sys.dim+1:end);
    Y = Eo(1:sys.dim,sys.dim+1:end);

else % MATLAB lyap solver

    P = lyap(Jr-Rr,(Jr-Rr)',Gr*Gr');
    Q = lyap((Jr-Rr)',Jr-Rr,Gr*Gr');
    X = lyap(A,(Jr-Rr)',G*Gr');
    Y = lyap(A',(Jr-Rr),-G*Gr');

end

%% Compute Gradient
grad.Jr = 2*((Q*P+Y'*X)-(Q*P+Y'*X).')/2;
grad.Rr = -2*Rr*(((Q*P+Y'*X)+(Q*P+Y'*X).')/2)*Rr;
grad.Gr = 2*((P+Q)*Gr + (Y-X)'*G);

end

function [hess] = hess(x,u,sys,redOrder,Opts)

Jr = x.Jr; Rr = x.Rr; Gr = x.Gr;
A = (sys.J-sys.R)*sys.Q; G = sys.G;
Jr_d = u.Jr; Rr_d = u.Rr; Gr_d = u.Gr;

if strcmp(Opts.lyapSolver,'mess')

    % Error system
    Ae = [A zeros(sys.dim,redOrder);zeros(redOrder,sys.dim) Jr-Rr];
    Be = [G;Gr];
    Ce = [G' -Gr'];

    % Solve Lyapunov equations of error system
    Z = mess_lyap(Ae, Be); Ec = Z*Z';
    Z = mess_lyap(Ae', Ce'); Eo = Z*Z';
    P = Ec(sys.dim+1:end,sys.dim+1:end);
    Q = Eo(sys.dim+1:end,sys.dim+1:end);
    X = Ec(1:sys.dim,sys.dim+1:end);
    Y = Eo(1:sys.dim,sys.dim+1:end);

    % Use workaround for coupled sparse-dense equations with
    % mess_lyap (Note: mess_lyap currently does not support solving
    % Sylvester equations)
    Be = [[X,G];[Jr_d-Rr_d,Gr_d]];
    Ce = [[Y';-G'],[(Jr_d-Rr_d);Gr_d']];
    Z = mess_lyap(Ae,Be);   Ec = Z*Z'; X_d =Ec(1:sys.dim,sys.dim+1:end);
    Z = mess_lyap(Ae',Ce'); Eo = Z*Z'; Y_d =Eo(1:sys.dim,sys.dim+1:end);
    % Use lyap for small, dense Sylvester equations
    P_d = lyap(Jr-Rr,(Jr-Rr)',(Jr_d-Rr_d)*P+P*(Jr_d-Rr_d)'+Gr_d*Gr'+Gr*Gr_d');
    Q_d = lyap((Jr-Rr)',Jr-Rr,(Jr_d-Rr_d)'*Q + Q*(Jr_d-Rr_d) + Gr_d*Gr'+Gr*Gr_d');

else
    P = lyap(Jr-Rr,(Jr-Rr)',Gr*Gr');
    Q = lyap((Jr-Rr)',Jr-Rr,Gr*Gr');
    X = lyap(A,(Jr-Rr)',G*Gr');
    Y = lyap(A',(Jr-Rr),-G*Gr');

    P_d = lyap(Jr-Rr,(Jr-Rr)',(Jr_d-Rr_d)*P+P*(Jr_d-Rr_d)'+Gr_d*Gr'+Gr*Gr_d');
    Q_d = lyap((Jr-Rr)',Jr-Rr,(Jr_d-Rr_d)'*Q + Q*(Jr_d-Rr_d) + Gr_d*Gr'+Gr*Gr_d');
    X_d = lyap(A,(Jr-Rr)',X*(Jr_d-Rr_d)'+G*Gr_d');
    Y_d = lyap(A',(Jr-Rr),Y*(Jr_d-Rr_d) - G*Gr_d');
end

% Riemannian Hessian
sym_df1S = ((Q*P+Y'*X)+(Q*P+Y'*X).')/2;
sk_hf1S =  ((Q_d*P+Q*P_d+Y_d'*X+Y'*X_d)-(Q_d*P+Q*P_d+Y_d'*X+Y'*X_d).')/2;
sym_hf1S = ((Q_d*P+Q*P_d+Y_d'*X+Y'*X_d)+(Q_d*P+Q*P_d+Y_d'*X+Y'*X_d).')/2;
sym_hf2S = ((Rr_d*sym_df1S*Rr)+(Rr_d*sym_df1S*Rr).')/2;

hess.Jr = 2*sk_hf1S;
hess.Rr = -2*Rr*sym_hf1S*Rr-2*sym_hf2S;
hess.Gr = 2*(P_d*Gr+P*Gr_d-X_d'*G+Q_d*Gr+Q*Gr_d+Y_d'*G);

end