function [sysr, specFac, redSpecFac] = sfmor(fctHandle, sys, varargin)
% SFMOR - Obtain a passive reduced order model by applying an sss-MOR method
%         to the spectral factor of system sys and computing a reduced
%         model for the associated reduced spectral factor
%
% Syntax:
%   sysr = SFMOR(fctHandle, sys)
%   sysr = SFMOR(fctHandle, sys, Opts)
%
% Description:
%       sysr = sfmor(fctHandle, sys) computes a reduced order model by
%       applying model reduction to the spectral
%       factor [1] of sys and subsequently reversing the construction of
%       the spectral factor.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - fctHandle:    function handle to the MOR function that is applied to
%                       the spectral factor
%       - sys:          sss (phs) object, containing passive LTI (port-Hamiltonian) system
%
%       *Optional Input Arguments:*
%       - additional inputs to @fctHandle except options
%       - Opts:  structure with execution parameters
%           - .makePH:                  Define if sysr should be pH or not
%                                       [{false} / true]
%           - .checkPassivity:          Check passivity of sys.
%                                       [{false} / true]
%           - .backtrafo:               Method to transform the spectral factor back.
%                                       [{popov} / prl]
%           - .tol:                     Tolerance if D has not full rank, D = tol*eye(m)
%                                       [{1e-12} / positive double]
%           - .lyap:                    Method to solve (reduced) Lyapunov equation
%                                       [{'lyap'} / 'mmess']
%           - .are:                     Method to solve (large-scale) ARE equation
%                                       [{'auto'} / 'mmess' / 'icare']
%           - .(func2str(fctHandle)):   Opts struct forwarded to MOR method in fcthandle
%           - .ss2phs.*:                Options that will be passed to ss2phs. Please refer to 
%                                       documentation of the respective algorithm.
%           - .phs.*:                   Options that will be passed to phs. Please refer to 
%                                       documentation of the respective algorithm.
%           - .samPassive.*:            Options that will be passed to samPassive.
%                                       Please refer to documentation of the respective algorithm.
%           - .areOpts.*:               Options that will be passed to mess_lrri.
%                                       Please refer to documentation of the respective algorithm.
%           - .lyapOpts.*:              Options that will be passed to mess_lradi.
%                                       Please refer to documentation of the respective algorithm.
%
% Output Arguments:
%       - sysr:         ssRed (phsRed) object, containing reduced LTI (port-Hamiltonian) system
%       - specFac:      sss object, containing spectral factor of sys
%       - redSpecFac:   ssRed object, containing spectral factor of sysr
%
% See Also:
%       phs, ss2phs, samPassive, demo_sfmor
%
% References:
%       [1] T. Breiten and B. Unger. “Passivity preserving model
%           reduction via spectral factorization." In: Automatica 142
%           (2022), p. 110368.
%       [2] B. C. Moore. Principal component analysis in linear systems: 
%           Controllability, observability, and model reduction. 
%           IEEE Transactions on Automatic Control, 26(1):17–32, 1981.
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
% Authors:      Maximilian Bonauer, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%               If you use the software package 'M-M.E.S.S.' [3], please refer
%               to its specific license.
%-----------------------------------------------------------------------

%% Input Parsing
[fctHandle,sys,D,dz,varargin,Opts] = parseInputs(fctHandle,sys,varargin{:});

%% Solve algebraic Riccati equation
if strcmp(Opts.are,'mmess')
    % Setup eqn struct for M.-M.E.S.S.
    R       = chol(D+D');
    Rinv    = R\eye(size(D));

    % eqn = struct('A_', sys.A, 'B1', sys.B*Rinv, 'C1', Rinv'*sys.C, 'U', -sys.B*Rinv,...
    %     'V', (Rinv'*sys.C)', 'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
    %     'haveUV', 1, 'type', 'T');
    eqn = struct('A_', sys.A -sys.B*Rinv*(Rinv'*sys.C), 'B1', sys.B*Rinv, 'C1', Rinv'*sys.C,...
        'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
        'haveUV', 0, 'type', 'T');

    % Operator manager for M.-M.E.S.S.
    oper = operatormanager(Opts.areOpts.oper);

    solX = mess_lrri(eqn, Opts.areOpts, oper);
    X = solX.Z*solX.Z';

    res = sys.A'*X+X*sys.A+(X*sys.B-sys.C')/(D+D')*(X*sys.B-sys.C')';
    err = max(abs(res(:)));

    if err > 1e-8
        warning('Computational issues with solution of algebraic Riccati equation.\n');
        warning('Error of algebraic Riccati equation is %.8f.\n',err);
    end
else
    [X,~,~,info] = icare(-full(sys.A),-full(sys.B),0,full(D+D'),full(sys.C'),'anti');
    err = norm(sys.A'*X+X*sys.A+(X*sys.B-sys.C')/(D+D')*(X*sys.B-sys.C')')/norm(X);
    if (info.Report>1)
        warning('Computational issues with solution of algebraic Riccati equation.\n');
        warning('Error of algebraic Riccati equation is %.8f.\n',err);
    end
end

%% Compute spectral factor
L = (D+D')^(-0.5)*(sys.C - sys.B'*X);
M = sqrtm(full(D+D'));

% LTL = -X*sys.A-sys.A'*X;
% L = lrcf(LTL,1e-12);
% L = L(1:size(sys.B,2),:);
% M = chol(D+D');

specFac = sss(sys.A, sys.B, L, M);

%% Model reduction of spectral factor
redSpecFac = fctHandle(specFac, varargin{:}, Opts.(func2str(fctHandle)));

% Make sure redSpecFac is explicit
redSpecFac.A = redSpecFac.E\redSpecFac.A;
redSpecFac.B = redSpecFac.E\redSpecFac.B;
redSpecFac.E = eye(size(redSpecFac.A));

% Check stability
eigRSF = eig(full(redSpecFac.A));
if ~(max(eigRSF) < 0)
    error('Reduced spectral factor is unstable.')
end

%% Solve Lyapunov equation
if strcmp(Opts.backtrafo,'prl')
    if strcmp(Opts.lyap,'mmess')
        % Setup eqn struct for M.-M.E.S.S.
        lyapEqn=struct('A_', redSpecFac.A, 'B', redSpecFac.B, 'C', redSpecFac.C, 'haveE', 0, 'type', 'T');

        % Operator manager for M.-M.E.S.S.
        oper = operatormanager(Opts.lyapOpts.oper);

        % Get adi shifts
        Opts.lyapOpts.shifts.p = mess_para(lyapEqn, Opts.lyapOpts, oper);

        % Low rank adi
        solXr = mess_lradi(lyapEqn, Opts.lyapOpts, oper);
        Xr = solXr.Z*solXr.Z';

    else % matlab solver lyap
        Xr = lyap(redSpecFac.A',redSpecFac.C'*redSpecFac.C);
    end
end

%% Construct reduced order model
if strcmp(Opts.backtrafo,'prl')
    Cr = redSpecFac.B'*Xr + redSpecFac.D'*redSpecFac.C;
    Dr = 0.5*(redSpecFac.D'*redSpecFac.D) + 0.5*(D-D');

    % Set feedthrough to zero if artificial D was introduced at the start
    if dz
        Dr = zeros(size(Dr));
    end

    % Reduced order model
    sysrSS = ssRed(redSpecFac.A, redSpecFac.B, Cr, Dr, eye(size(redSpecFac.A)), 'sfmor', Opts);

else

    % Obtain original system via Popov function
    SpF = ss(redSpecFac.A, redSpecFac.B, redSpecFac.C, redSpecFac.D);
    PP = SpF'*SpF;
    Gs = stabsep(PP);
    Gs.D = Gs.D/2;

    Ar = Gs.A;
    Br = Gs.B;
    Cr = Gs.C;
    Dr = Gs.D;

    % Set feedthrough to zero if artificial D was introduced at the start
    if dz
        Dr = zeros(size(Dr));
    end

    sysrSS = ssRed(Ar,Br,Cr,Dr,eye(size(Ar)),'sfmor',Opts);
end

% Make reduced pH system
if Opts.makePH
    if strcmp(Opts.backtrafo,'prl')
        sysr = ss2phs(sysrSS, Xr, Opts.ss2phs);
    else
        sysr = ss2phs(sysrSS, Opts.ss2phs);
    end
    sysr.method = @sfmor;
    sysr.parameters = Opts;
    sysr.parameters.sysSS = sysrSS;
else
    sysr = sysrSS;
end

end

%% Auxilliary functions
function [fctHandle,sys,D,dz,varargin,Opts] = parseInputs(fctHandle,sys,varargin)
if nargin < 2
    error('Not enough input arguments.')
end

if ~isa(fctHandle,'function_handle')
    error('First argument is not a function_handle.');
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end

% Set default options
OptsAdmissible.(func2str(fctHandle)) = struct();
OptsAdmissible.makePH = {false,true};
OptsAdmissible.backtrafo = {'popov','prl'};
OptsAdmissible.tol = 1e-12;
OptsAdmissible.lyap = {'lyap','mmess'};                 % Use lyap for Lyapunov equation
OptsAdmissible.are = {'auto','mmess','icare'};          % Toolbox to solve ARE
OptsAdmissible.checkPassivity = {false,true};
OptsAdmissible.phs = struct();                          % Make sure options struct for phs exists
OptsAdmissible.ss2phs = struct();                       % Make sure options struct for ss2phs exists
OptsAdmissible.samPassive.plot = false;                 % Do not plot Popov function

% Set mmess ARE options
OptsAdmissible.areOpts.norm = {'fro',2};
OptsAdmissible.areOpts.ri.riccati_solver = 'radi';
OptsAdmissible.areOpts.ri.lqg_solver = 'radi';
OptsAdmissible.areOpts.ri.maxiter = 500;
OptsAdmissible.areOpts.ri.res_tol = 1e-15;
OptsAdmissible.areOpts.ri.rel_diff_tol = 1e-15;
OptsAdmissible.areOpts.ri.trunc_tol = 1e-15;
OptsAdmissible.areOpts.ri.compres_tol = 1e-15;
OptsAdmissible.areOpts.ri.info = 0;
OptsAdmissible.areOpts.ri.Z0 = [];
OptsAdmissible.areOpts.ri.store_lqg = 0;
OptsAdmissible.areOpts.ri.store_solfac = 0;
OptsAdmissible.areOpts.ri.trunc_info = 0;
OptsAdmissible.areOpts.shifts.num_desired = 10;
OptsAdmissible.areOpts.shifts.method = 'projection';
OptsAdmissible.areOpts.adi.maxiter = 500;
OptsAdmissible.areOpts.adi.res_tol = 1e-15;
OptsAdmissible.areOpts.adi.rel_diff_tol = 1e-15;
OptsAdmissible.areOpts.adi.info = 0;
OptsAdmissible.areOpts.adi.compute_sol_fac = 0;
OptsAdmissible.areOpts.adi.accumulateK = 0;
OptsAdmissible.areOpts.adi.accumulateDeltaK = 0;
OptsAdmissible.areOpts.radi.maxiter = 500;
OptsAdmissible.areOpts.radi.res_tol = 1e-15;
OptsAdmissible.areOpts.radi.rel_diff_tol = 1e-15;
OptsAdmissible.areOpts.radi.info = 0;
OptsAdmissible.areOpts.oper = 'default';

% Set mmess lyap options
OptsAdmissible.lyapOpts.adi.maxiter = 500;
OptsAdmissible.lyapOpts.adi.res_tol = 0;
OptsAdmissible.lyapOpts.adi.rel_diff_tol = 0;
OptsAdmissible.lyapOpts.adi.info = 0;
OptsAdmissible.lyapOpts.shifts.num_desired = 25;
OptsAdmissible.lyapOpts.shifts.num_Ritz = 50;
OptsAdmissible.lyapOpts.shifts.num_hRitz = 25;
OptsAdmissible.lyapOpts.shifts.b0 = ones(size(sys.A,1),1);
OptsAdmissible.lyapOpts.shifts.method = 'projection';
OptsAdmissible.lyapOpts.shifts.info = 0;
OptsAdmissible.lyapOpts.norm = {'fro',2};
OptsAdmissible.lyapOpts.oper = 'default';

Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

% Set values for ss2phs
Opts.ss2phs.phs = Opts.phs;
Opts.ss2phs.method = @sfmor;

% Check for pH system
isPHS = 0;
if isa(sys,'phs') || isa(sys,'phsRed')
    sys = phs2sss(sys);
    isPHS = 1;
end

% Check for passive system
if ~isPHS && Opts.checkPassivity
    passive = samPassive(sys, Opts.samPassive);
else
    passive = true;
end

if ~passive
    error('Spectral factor MOR is only applicable to passive systems.')
end

% Check for D+D' > 0 and correct if necessary
ev = eig(sys.D+sys.D');
dz = false;
if all(ev>0)
    % Positive Definite
    D = sys.D;
else
    D = sys.D + Opts.tol*eye(size(sys.D));
    dz = true;
end

% Check if sys is DAE
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('sfmor currently only supports ODE systems.')
end

% Check if E exists
if ~isempty(sys.E)
    % Make sure sys is explicit
    sys.A = sys.E\sys.A;
    sys.B = sys.E\sys.B;
    sys.E = eye(size(sys.A));
else
    sys.E = eye(size(sys.A));
end

% Select ARE solver if 'auto'
if strcmp(Opts.are,'auto')
    if sys.n > 1000
        Opts.are = 'mmess';
    else
        Opts.are = 'icare';
    end
end

% Check third-party software
if strcmp(Opts.lyap,'mmess') || strcmp(Opts.are,'mmess')
    thirdPartyCheck('M.E.S.S.');
end

end
