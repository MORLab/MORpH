function [redSys] = prOpt(sys,redOrder,varargin)
% PROPT - Optimizes the matrix entries of a reduced-order model
% in an H2-optimal way using the pole-residue-framework
%
% Syntax:
%   redSys = PROPT(sys,redOrder)
%   redSys = PROPT(sys,redOrder,redSys0)
%   redSys = PROPT(sys,redOrder,Opts)
%   redSys = PROPT(sys,redOrder,redSys0,Opts)
%
% Description:
%       redSys = prOpt(sys, redOrder) returns a reduced PH system of order
%       redOrder.
%
%       Argument Opts can be used for further adjustments.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object (large-scale model)
%       - redOrder:	desired order of the proper part of the reduced-order model
%                   If an initial model redSys0 is provided, its proper
%                   part should have the same dimension as redOrder
%
%       *Optional Input Arguments:*
%       - redSys0: Initial phs object used as a starting point of the
%                  optimization
%       - Opts:  structure with execution parameters
%           - .optParams:       Structure that specifies which matrices shall be optimized
%                               Default: struct('Er',false,'Jr',true,'Rr',true,'Qr',true,'Gr',true,'Pr',true)
%           - .eigDerivative:   Way of computing the derivatives of the
%                               eigenvalue problem. 
%                               [{'magnus'} / 'rudisill' / 'murthy']
%           - .GspMethod:       Determines whether the strictly proper part of G(s) is sampled directly or
%                               indirectly via Gsp(s) = G(s)-Pol(s)
%                               [{'indirect'} / 'direct']
%           - .maxIter:         Defines the maximum amount of iterations
%                               [{1000} / positive integer]
%           - .tol:             Optimality tolerance (varies for each
%                               solver)
%                               [{1e-7} / positive double]
%           - .solver:       	Optimization solver
%                               [{'Manopt'} / 'GRANSO' / 'MATLAB']
%           - .decomposePHDAE.*:Other options that will be passed on to the function;
%                               Please refer to its documentation (doc decomposePHDAE).
%           - .composePHDAE.*:  Other options that will be passed on to the function;
%                               Please refer to its documentation (doc composePHDAE).
% Output Arguments:
%       - redSys: 	phsRed object (reduced-order model)
%
% See Also:
%       dPR_dx, lyapOpt, composePHDAE, decomposePHDAE
%
% References:
%       [1] T. Moser and B. Lohmann. “A New Riemannian Framework for Efficient
%           H2-Optimal Model Reduction of Port-Hamiltonian Systems." 
%           In: Proceedings of 59th IEEE Conference on Decisison and Control (CDC). 
%           Jeju Island, Republic of Korea, 2020, pp. 5043–5049.
%       [2] T. Moser et al. Structure-Preserving Model Order Reduction for Index 
%           Two Port-Hamiltonian Descriptor Systems. 
%           arXiv Preprint arXiv:2206.03942. 2022. url: https://arxiv.org/abs/2206.03942
%       [3] P. Schwerdtner et al. Structure-Preserving Model Order
%           Reduction for Index One Port-Hamiltonian Descriptor Systems. 
%           arXiv Preprint arXiv:2206.01608. 2022. url: https://arxiv.org/abs/2206.01608.
%       [4] N. Boumal et al. "Manopt, a Matlab Toolbox for Optimization on Manifolds." 
%           In: Journal of Machine Learning Research 15.42 (2014), pp. 1455–1459. 
%           url: https://www.manopt.org.
%       [5] F. E. Curtis, T. Mitchell, and M. L. Overton. "A BFGSSQP method for 
%           nonsmooth, nonconvex, constrained optimization and its evaluation 
%           using relative minimization profiles.” In: Optimization Methods 
%           and Software 32.1 (2017), pp. 148–181.
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
%               If you use the software packages 'Manopt' [4] or 'GRANSO' [5], 
%               please refer to their specific license.
%-----------------------------------------------------------------------

%% Input parsing
[redSys0, Opts] = parseInputs(sys, redOrder, varargin{:});

% Get dimensions
nInp = size(sys.G,2);

%% Decompose original transfer function
if sys.isDAE
    % Decomposition in proper and improper part
    [sysP, Pol_1] = decomposePHDAE(sys,Opts.decomposePHDAE);
    if any(eig(Pol_1)<0)
        if min(eig(Pol_1)>-1e-12)
            warning('MORpH:prOpt:decomposePHDAE',"Enforcement of positive semidefinite polynomial part may lead to unbounded H2 errors.")
            Pol_1 = makePSD(Pol_1);
        else
            error('MORpH:prOpt:decomposePHDAE',"Computed linear polynomial part of original system is negative definite");
        end
    end
else
    Pol_1 = zeros(size(sys.G,2));
    sysP = sys;
end

% Compute data to evaluate the original transfer function
FOM_eval = struct();
if strcmp(Opts.GspMethod,'direct')
    sysSP = phs2sss(sysP);
    sysSP.D = zeros(size(sysSP.D));
    FOM_eval.sysSP = sysSP;
else
    FOM_eval.sys = sys;
    FOM_eval.Pol_0 = sysP.S+sysP.N;
    FOM_eval.Pol_1 = Pol_1;
end

% Compute feedthrough matrices Sr and Nr
Nr = full(sysP.N);
Sr = full(sysP.S);

if isequal(Sr,zeros(nInp,nInp))
    % If Sr is zero, Pr has to be as well
    Opts.optParams.Pr = false;
    Ls = zeros(nInp,nInp);
else
    [Us,Ds] = eig(full(Sr));
    Ls = Us*sqrt(Ds);
end

%% Initialization
% Create the problem structure
params = fieldnames(Opts.optParams);
dims = struct('nx',redOrder,'nu',nInp,'nE',0,'nJ',0,'nR',0,'nQ',0,'nG',0,'nP',0,'nL',0);
for im = 1:numel(params)
    if Opts.optParams.(params{im})
        switch params{im}
            case 'Er'
                dims.nE = dims.nx*(dims.nx+1)/2;
            case 'Jr'
                dims.nJ = dims.nx*(dims.nx-1)/2;
            case 'Rr'
                dims.nR = dims.nx*(dims.nx+1)/2;
            case 'Qr'
                dims.nQ = dims.nx*(dims.nx+1)/2;
            case 'Gr'
                dims.nG = dims.nx*dims.nu;
            case 'Pr'
                dims.nP = dims.nx*dims.nu;
        end
    end
end
dims.nvar = dims.nE+dims.nJ+dims.nR+dims.nQ+dims.nG+dims.nP;

% Create initial parameter theta0 of optimization parameters and struct
% 'fix' with fixed matrices

if isa(redSys0,'phs')
    [theta0, fix] = parseInitialROM(redSys0,dims,sysP,Pol_1,Opts);

else
    % random initialization
    if dims.nE, theta0E = utv(eye(redOrder));                 else, theta0E = []; fix.Er = eye(dims.nx); end
    if dims.nJ, theta0J = sutv(triu(-1+2*rand(redOrder),1));  else, theta0J = []; fix.Jr = zeros(dims.nx); end
    if dims.nR, theta0R = utv(eye(redOrder));                 else, theta0R = []; fix.Lr = zeros(dims.nx); end
    if dims.nQ, theta0Q = utv(eye(redOrder));                 else, theta0Q = []; fix.Qr = eye(dims.nx); end
    if dims.nG, theta0G = ftv(-1+2*rand(redOrder,nInp));      else, theta0G = []; fix.Gr = zeros(dims.nx,dims.nu); end
    if dims.nP, theta0P = ftv(zeros(redOrder,nInp));          else, theta0P = []; fix.Lp = zeros(dims.nx,dims.nu); end

    theta0 = [theta0E;theta0J;theta0R;theta0Q;theta0G;theta0P];

end

% Fix feedthrough matrices
fix.Ls = Ls;
fix.Nr = Nr;

%% Optimization
warning('off','MORpH:phs:noInputValidation')
if Opts.maxIter > 0
    switch Opts.solver
        case 'Manopt'
            problem.M = euclideanfactory(dims.nvar,1);
            problem.costgrad = @(theta) objectiveFunction(theta,fix,dims,FOM_eval,Opts);
            if Opts.manopt.checkGradient
                checkgradient(problem);
            end
            if Opts.manopt.checkHessian
                checkhessian(problem);
            end
            options.tolgradnorm = Opts.tol;
            options.maxiter = Opts.maxIter;
            options.minstepsize = 1e-20;
            [theta_opt, ~, info] = trustregions(problem,theta0,options);
        case 'GRANSO'
            % Create the objective function for GRANSO
            objFun = @(theta)objectiveFunction(theta,fix,dims,FOM_eval,Opts);
            % Start optimization with GRANSO
            options.x0 = theta0;
            options.maxit = Opts.maxIter;
            options.opt_tol = Opts.tol;
            solution = granso(dims.nvar, objFun, [], [], options);
            theta_opt = solution.best.x;
            info = [];
        case 'MATLAB'
            % Create the objective funciton for fminunc
            objFun = @(theta)objectiveFunction(theta,fix,dims,FOM_eval,Opts);
            % Start optimization with fminunc
            options = optimoptions(@fminunc,'Algorithm','trust-region',...
                'SpecifyObjectiveGradient',true,'Display','iter', ...
                'OptimalityTolerance',Opts.tol,'MaxIterations',Opts.maxIter);
            theta_opt = fminunc(objFun,theta0,options);
            info = [];
        otherwise
    end

else
    warning('MORpH:prOpt:wrongInput',['The desired maximum number of iterations is smaller than 1.\n',...
        'Exit with the initial model.']);
    theta_opt = theta0;
    info = [];
end
warning('off','MORpH:phs:noInputValidation');

%% Create redSys object
redSysP = create_ROM(theta_opt, fix, dims, Opts);
% Attach feedthrough matrices
redSysP.N = Nr;
redSysP.S = Sr;

% Attach polynomial part (if necessary)
redSys = composePHDAE(redSysP, Pol_1, Opts.composePHDAE);

% Input validation and storage of additional information
phs.inputValidation(redSys);
redSys.method = @prOpt;
redSys.parameters = Opts;
redSys.info.theta = theta_opt;
redSys.info.init = theta0;
redSys.info.convergence = info;

end

%% Auxiliary functions

function [redSys0,Opts] = parseInputs(sys, redOrder, varargin)

% check number of inputs provided
narginchk(2,4);

% Check original system
if ~isa(sys,'phs')
    error('MORpH:prOpt:wrongInput', ...
        ['Original model should be a pHDAE of the following form:\n',...
        'E*dx/dt = (J-R)*Q*x(t) + (G-P)*u(t)\n',...
        'y(t) = (G+P)^T*Q*x(t)+(S+N)*u(t)']);
end

% Check redOrder
if redOrder > sys.dim
    error('MORpH:prOpt:wrongInput', 'Reduced Order should be smaller than original system order')
end

% Parse optional arguments
if nargin==2
    Opts = struct();
    redSys0 = false;
elseif nargin==3
    if isstruct(varargin{1})
        Opts = varargin{1};
        redSys0 = false;
    elseif isa(varargin{1},'phs')
        redSys0 = varargin{1};
        Opts = struct();
    else
        error('MORpH:prOpt:wrongInput', 'Wrong input combination');
    end
else
    redSys0 = varargin{1};
    Opts = varargin{2};
end

% Check if opts are admissible
optsAdmissible.optParams = struct();
optsAdmissible.eigDerivative = {'magnus','murthy','rudisill'};
optsAdmissible.GspMethod = {'indirect','direct'};
optsAdmissible.maxIter = 1000;
optsAdmissible.tol = 1e-7;
optsAdmissible.solver = {'Manopt','GRANSO','MATLAB'};
optsAdmissible.manopt.checkGradient = {false,true};
optsAdmissible.manopt.checkHessian = {false,true};
optsAdmissible.phs.inputValidation = false; % ROM is checked at the end
optsAdmissible.composePHDAE.type = 'phsRed';
optsAdmissible.decomposePHDAE = struct();

Opts = phsMOR_parseOpts(Opts,optsAdmissible);

% Set non-specified optimization parameters to true
optParamsNew = struct('Er',false,'Jr',true,'Rr',true,'Qr',true,'Gr',true,'Pr',true);
optParams = fieldnames(optParamsNew);
for i = 1:numel(optParams)
    if isfield(Opts.optParams,(optParams{i})) && isa(Opts.optParams.(optParams{i}),'logical')
        optParamsNew.(optParams{i}) = Opts.optParams.(optParams{i});
    end
end
Opts.optParams = optParamsNew;
assert(~(Opts.optParams.Er && Opts.optParams.Qr), "It is not meaningful to optimize Er and Qr at the same time " + ...
    "because they equally contribute to the Hamiltonian. Please optimize either of them.");

% Check for 3rd-party software
switch Opts.solver
    case 'Manopt'
        thirdPartyCheck('Manopt');
    case 'GRANSO'
        thirdPartyCheck('GRANSO');
end

end

function [theta0, fix] = parseInitialROM(redSys0, dims, sysP, Pol_1, Opts)

if redSys0.isDAE
    % Decomposition in proper and improper part
    [redSys0P, redPol_1] = decomposePHDAE(redSys0,Opts.decomposePHDAE);
    if any(eig(Pol_1)<0)
        if min(eig(Pol_1))>-1e-12
            warning('MORpH:prOpt:decomposePHDAE',"Enforcement of positive semidefinite polynomial part may lead to unbounded H2 errors.")
            Pol_1 = makePSD(Pol_1,1e-12);
        else
            error('MORpH:prOpt:decomposePHDAE',"Computed linear polynomial part of initial system is negative definite");
        end
    end
else
    redPol_1 = zeros(size(redSys0.G,2));
    redSys0P = redSys0;
end

% Check initial reduced model
if redSys0P.dim ~= dims.nx
    error('MORpH:prOpt:wrongInput',['The specified reduced order and the dimension ' ...
        'of the proper part of the provided initial model do not match']);
end
if isa(redSys0,'phs')
    if size(redSys0P.G,2) ~= size(sysP.G,2)
        error('MORpH:prOpt:wrongInput', ...
            'Initial model should have the same number of inputs as the original model.');
    end
end

% Ensure that Q'*E = E'*Q > 0 is always true
if Opts.optParams.Qr && ~isequal(redSys0P.E,eye(size(redSys0P.E)))
    redSys0P = makeExplicit(redSys0P);
elseif Opts.optParams.Er && ~isequal(redSys0P.Q,eye(size(redSys0P.Q)))
    redSys0P = scaling(redSys0P);
end

% Check if sys and redSys0 have different polynomial parts
if norm(redPol_1 - Pol_1) > 1e-12
    warning("Original and initial model have different improper parts. " + ...
        "We proceed with the improper part of the original model.");
end
if norm(full(sysP.N-redSys0P.N)) > 1e-12
    warning("Original and initial model have different feedthrough matrices N. " + ...
        "We proceed with the matrix N of the original model.");
    redSys0P.N = sysP.N;
end
if norm(full(sysP.S-redSys0P.S)) > 1e-12
    redSys0P.S = sysP.S;
    redSys0P.P = zeros(size(redSys0P.P));
    warning("Original and initial model have different feedthrough matrices S. " + ...
        "We proceed with the matrix S of the original model and initialize with P==0.");
end

% Make R0 or W0 positive definite if not already the case
nInp = size(redSys0.S,2);
if isequal(redSys0.S,zeros(nInp))
    % Make Rr positive definite (if necessary)
    Rpd = makePSD(full(redSys0P.R),1e-12);
    % Compute Cholesky factors
    Lr = chol(full(Rpd),'lower');
    Lp = zeros(size(redSys0P.P));
else
    W = full([redSys0P.S,redSys0P.P';redSys0P.P,redSys0P.R]);
    % Make W positive definite (if necessary)
    Wpd = makePSD(W,1e-12);
    % Compute Cholesky factors
    L = chol(Wpd,'lower');
    Lr = L(nInp+1:end,nInp+1:end);
    Lp = L(nInp+1:end,1:nInp);
end

% Compute theta0 and fixed matrices
if dims.nE, theta0E = utv(chol(full(redSys0P.E)));        else, theta0E = []; fix.Er = redSys0P.E; end
if dims.nJ, theta0J = -sutv(triu(full(redSys0P.J),1));    else, theta0J = []; fix.Jr = redSys0P.J; end
if dims.nR, theta0R = utv(Lr.');                          else, theta0R = []; fix.Lr = Lr; end
if dims.nQ, theta0Q = utv(chol(full(redSys0P.Q)));        else, theta0Q = []; fix.Qr = redSys0P.Q; end
if dims.nG, theta0G = ftv(redSys0P.G);                    else, theta0G = []; fix.Gr = redSys0P.G; end
if dims.nP, theta0P = ftv(Lp);                            else, theta0P = []; fix.Lp = Lp; end
theta0 = [theta0E;theta0J;theta0R;theta0Q;theta0G;theta0P];

end

function [f, grad] = objectiveFunction(theta, fix, dims, FOM_eval, Opts)

% Create full ROM
redSys = create_ROM(theta,fix,dims,Opts);
redOrder = redSys.dim;
nInp = size(redSys.G,2);

% Strictly proper part of Gr(s)
Ersp = redSys.E;
Arsp = (redSys.J-redSys.R)*redSys.Q;
Brsp = redSys.G-redSys.P;
Crsp = (redSys.G+redSys.P)'*redSys.Q;

% Compute eigenvector for each lambdaR
[Xr,Dr,Xt] = eig(Arsp,Ersp);
% Normalize eigenvectors
%         Xr = Xr./vecnorm(Xr);

% Sort eigenvectors and -values according to lambdaR
[~,ind]=sort(diag(Dr),'descend');
Dr = Dr(ind,ind);
Xr = Xr(:,ind);
Xt = Xt(:,ind);
lambdaR = diag(Dr);
Yr = Xr\eye(redOrder);

switch Opts.eigDerivative
    case 'murthy'
        % Determine row m of max element
        [~,m_max] = max(abs(Yr).'.*abs(Xr),[],1);
        Gamma = diag(1./diag(Xr(m_max,:)));
        % Transformation of EVP
        Xr = Xr*Gamma;
    otherwise
end

% Matrices for modal coordinate residuals
Lt = Crsp*Xr;
Rt = Xr\(Ersp\Brsp);


%     % Check for duplicate eigenvalues
%     if length(lambdaR)-length(unique(lambdaR))>0
%         error('Duplicate eigenvalues');
%     end

% Evaluate Gsp, dGsp at -lambdaR
Gsp=zeros(nInp,nInp,length(lambdaR));
dGsp=zeros(nInp,nInp,length(lambdaR));

for j=1:redOrder
    [Gsp(:,:,j),dGsp(:,:,j)] = eval_Gsp(-lambdaR(j),FOM_eval,Opts.GspMethod);
end

%% Cost function

f = 0;
for k=1:redOrder
    f = f - 2*Lt(:,k).'*Gsp(:,:,k)*Rt(k,:).' + Lt(:,k).'*(Crsp*((-lambdaR(k)*Ersp-Arsp)\Brsp))*Rt(k,:).';
end
f = real(f);

%% Gradient

df_dLambdaR = zeros(1,redOrder);
df_dPhiR = zeros(1,2*redOrder*nInp);
j=1;
ptr=1;
while j<=redOrder

    df_dLambdaR(j) = -2*Lt(:,j).'*(-Crsp*(((-lambdaR(j)*eye(redOrder,redOrder)-(Ersp\Arsp))*(-lambdaR(j)*eye(redOrder,redOrder)-(Ersp\Arsp)))\(Ersp\Brsp)) - dGsp(:,:,j))*Rt(j,:).';
    xi_inv = Crsp*((-lambdaR(j)*eye(redOrder,redOrder)-(Ersp\Arsp))\(Ersp\Brsp)) - Gsp(:,:,j);
    df_dPhiR(ptr:(ptr+nInp-1)) = (2*xi_inv*Rt(j,:).').';
    df_dPhiR((ptr+nInp):(ptr+2*nInp-1)) = 2*Lt(:,j).'*xi_inv;

    if j<redOrder && lambdaR(j+1)==conj(lambdaR(j))
        df_dLambdaR(j+1) = conj(df_dLambdaR(j));
        df_dPhiR((ptr+2*nInp):(ptr+4*nInp-1)) = conj(df_dPhiR(ptr:(ptr+2*nInp-1)));
        j=j+2;
        ptr = ptr+4*nInp;
    else
        j=j+1;
        ptr = ptr+2*nInp;
    end
end

% Derivative vector (q = [phiR(r);lambdaR(r)]
df_dq = [df_dPhiR,df_dLambdaR];

% Derivative of PhiR, LambdaR w.r.t. theta
evp.Xr = Xr;
evp.Xl = Xt;
evp.Lambda = Dr;
[dPhiR_dx,dLambdaR_dx] = dPR_dx(redSys,theta,dims,evp,Opts.eigDerivative);

J_ = [dPhiR_dx;dLambdaR_dx];

grad = real(df_dq*J_)';
end

function [Gsp_si, dGsp_si] = eval_Gsp(s_i, FOM, method)
% Evaluates the strictly-proper transfer function Gsp(s) and its
% derivative dGsp(s)/ds at certain shifts s_i. This happens either
% directly with the system matrices of the strictly-proper part or
% indirectly using the fact that Gsp(s) = G(s)-Pol(s)

if strcmp(method,'direct')
    res1 = (s_i*FOM.sysSP.E-FOM.sysSP.A)\FOM.sysSP.B;
    Gsp_si = FOM.sysSP.C*res1;
    dGsp_si = -FOM.sysSP.C*((s_i*FOM.sysSP.E-FOM.sysSP.A)\(FOM.sysSP.E*res1));

else % Gsp(s) = G(s)-Pol(s)
    res1 = (s_i*FOM.sys.E-(FOM.sys.J-FOM.sys.R)*FOM.sys.Q)\(FOM.sys.G-FOM.sys.P);
    Gsp_si = (FOM.sys.G+FOM.sys.P)'*FOM.sys.Q*res1 + (FOM.sys.S+FOM.sys.N) - (FOM.Pol_0 + FOM.Pol_1*s_i);
    dGsp_si = -(FOM.sys.G+FOM.sys.P)'*FOM.sys.Q*((s_i*FOM.sys.E-(FOM.sys.J-FOM.sys.R)*FOM.sys.Q)\(FOM.sys.E*res1)) - FOM.Pol_1;
end
end

function [thetaE, thetaJ, thetaR, thetaQ, thetaG, thetaP] = unziptheta(theta, dims)
% Decompose parameter vector
thetaE = theta(1:dims.nE); theta = theta(dims.nE+1:end);
thetaJ = theta(1:dims.nJ); theta = theta(dims.nJ+1:end);
thetaR = theta(1:dims.nR); theta = theta(dims.nR+1:end);
thetaQ = theta(1:dims.nQ); theta = theta(dims.nQ+1:end);
thetaG = theta(1:dims.nG); theta = theta(dims.nG+1:end);
thetaP = theta(1:dims.nP);
end

function sys = create_ROM(theta, fix, dims, opts)
[thetaE, thetaJ, thetaR, thetaQ, thetaG, thetaP] = unziptheta(theta, dims);

% Construct system matrices
if dims.nE, E = vtu(thetaE)'*vtu(thetaE);      else, E = fix.Er;   end
if dims.nJ, J = vtsu(thetaJ)'-vtsu(thetaJ);    else, J = fix.Jr; end
if dims.nR, R = vtu(thetaR)'*vtu(thetaR);      else, R = fix.Lr*fix.Lr'; end
if dims.nQ, Q = vtu(thetaQ)'*vtu(thetaQ);      else, Q = fix.Qr;   end
if dims.nG, G = vtf(thetaG, dims.nx, dims.nu); else, G = fix.Gr; end
if dims.nP
    P = vtf(thetaP, dims.nx, dims.nu)*fix.Ls';
    R = R + vtf(thetaP, dims.nx, dims.nu)*vtf(thetaP, dims.nx, dims.nu)';
else
    P = fix.Lp*fix.Ls';
    R = R + fix.Lp*fix.Lp';
end
sys = phsRed(J, R, Q, G, E, P, fix.Ls*fix.Ls', fix.Nr, opts.phs);

end