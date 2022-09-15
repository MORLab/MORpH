function sysPH = ss2phs(sys, varargin)
% ss2phs - Returns pH system representation of a given passive system
%
% Syntax:
%   sysPH = ss2phs(sys)
%   sysPH = ss2phs(sys,X)
%   sysPH = ss2phs(sys,Opts)
%   sysPH = ss2phs(sys,X,Opts)
%
% Description:
%       This funciton uses the Kalman-Yakubovic-Popov Lemma to construct a
%       port-Hamiltonian system representation of a given passive
%       ss object.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  ss object to be brought in pH representation
%       *Optional Input Arguments:*
%       - X: solution to positive real Riccati equation
%       - Opts: structure with execution parameters
%           - .tol:             If D+D' is not positive definite, D is replaced
%                               by Opts.tol*eye(m).
%                               [{1e-12} / positive double]
%           - .representation:  Defines the pH system representation
%                               [{'implicit'} / 'standard' / 'scaled']
%           - .method:          Reduction method assigned to phsRed object.
%                               [{@ss2phs} / function handle]
%           - .prl:             Define how to solve the Positive Real lemma.
%                               [{'cvx'} / 'icare' / 'mmess' / 'yalmip']
%           - .samPassive.*:    Other options that will be passed on to the
%                               method samPassive;
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .phs.*:           Other options that will be passed on to the
%                               class 'phs';
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .mmess.*:         Other options that will be passed on to mess_lrri
%                               of the M-M.E.S.S. toolbox to solve the PRL ARE.
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .yalmip.solver:   Solver option that will be passed on to sdpsettings
%                               to solve the PRL LMI.
%                               Refer to documentation of the respective
%                               algorithm.
%                               [{'sedumi'} / string]
%
% Output Arguments:
%       - sysPH: 	phsRed object
%
% See Also:
%       samPassive, phs, mess_lrri, icare, cvx, yalmip
%
% References:
%       [1] J. Saak, M. Köhler, and P. Benner. M-M.E.S.S.-2.1 – The
%           Matrix Equations Sparse Solvers library.
%           URL: https://www.mpi-magdeburg.mpg.de/projects/mess
%       [2] M. Grant and S. Boyd. CVX: Matlab Software for Disciplined Convex 
%           Programming, version 2.1. Mar. 2014. url: http://cvxr.com/cvx.
%       [3] M. Grant and S. Boyd. “Graph implementations for
%           nonsmooth convex programs." In: Recent Advances in
%           Learning and Control. Ed. by V. Blondel, S. Boyd, and H.
%           Kimura. Lecture Notes in Control and Information Sciences.
%           Springer-Verlag Limited, 2008, pp. 95–110.
%       [4] J. Löfberg. “YALMIP : A Toolbox for Modeling and Optimization in MATLAB.” 
%           In: In Proceedings of the CACSD Conference. Taipei, Taiwan, 2004.
%       [5] C. Beattie, V. Mehrmann, and P. Van Dooren. Robust port-Hamiltonian 
%           representations of passive systems. Automatica J. IFAC, 100:182–186, 2019.
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
%               If you use the software packages 'M-M.E.S.S.' [1], 'CVX' [2,3] or 'YALMIP' [4], 
%               please refer to their specific license.
%-----------------------------------------------------------------------

%% Input Parsing
[sys,D,X,availableX,Opts] = parseInputs(sys,varargin{:});

%% Minimal systems
if ~availableX
    switch Opts.prl
        case 'mmess'
            % Setup eqn struct for M.-M.E.S.S.
            Z       = chol(D+D');
            Zinv    = Z\eye(size(D));

            %             eqn = struct('A_', sys.A, 'B1', sys.B*Zinv, 'C1', Zinv'*sys.C, 'U', -sys.B*Zinv,...
            %                 'V', (Zinv'*sys.C)', 'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
            %                 'haveUV', 1, 'type', 'T');
            eqn = struct('A_', sys.A-sys.B*Zinv*(Zinv'*sys.C), 'B1', sys.B*Zinv, 'C1', Zinv'*sys.C,...
                'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
                'haveUV', 0, 'type', 'T');

            % Operator manager for M.-M.E.S.S.
            oper = operatormanager(Opts.mmess.oper);

            % Minimal solution to Riccati equation
            sol = mess_lrri(eqn, Opts.mmess, oper);
            X = sol.Z*sol.Z';

        case 'icare'
            [X,~,~,info] = icare(full(sys.A),full(sys.B),0,full(-D-D'),full(-sys.C'));

            if info.Report == 1
                warning('Solution accuracy to KYP inequality is poor.')
            elseif info.Report > 1
                error('No accurate solution to KYP inequality found.')
            end

        case 'cvx'
            n = size(sys.A,1);
            m = size(sys.B,2);

            cvx_begin quiet sdp
            variable X(n,n);
            -[sys.A'*X+X'*sys.A, X*sys.B-sys.C'; sys.B'*X-sys.C, -D-D'] == semidefinite(n+m);
            X - 1e-8*eye(n) == semidefinite(n);
            cvx_end

            if strcmp(cvx_status,'Solved')
                % Nothing to do
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                warning('Solution to KYP inequality is inaccurate.')
            else
                error('No accurate solution to KYP inequality found.')
            end

        case 'yalmip'
            yalmip('clear')
            n = size(sys.A,1);

            opt = sdpsettings('solver',Opts.yalmip.solver,'verbose',0); % Optimization settings

            % Define variables
            X = sdpvar(n,n);

            % Define matrices
            mat = [sys.A'*X+X'*sys.A, X*sys.B-sys.C'; sys.B'*X-sys.C, -D-D'];

            % Define Constraints
            F = [];
            F = [F X>=1e-8*eye(n)];
            F = [F mat<=-1e-8*eye(size(mat))];

            % Optimization Problem
            info = optimize(F,[],opt); 
            
            X = value(X);
            
            if info.problem == 0
                % Nothing to do
            elseif info.problem == 1
                warning('YALMIP thinks problem is infeasible.')
            else
                warning('Solution accuracy to KYP inequality may be inaccurate.')
            end

            if any(isnan(X))
                error('No accurate solution to KYP inequality found.')
            end
    % End of switch Opts.prl        
    end

    % End of if ~availableX
end

% Construct system matrices
switch Opts.representation
    case 'standard'
        Q = X;
        E = eye(size(Q));
        J = 0.5*(sys.A/X - X\sys.A');
        R = -0.5*(sys.A/X + X\sys.A');
        G = 0.5*(X\sys.C' + sys.B);
        P = 0.5*(X\sys.C' - sys.B);

    case 'scaled'
        T = chol(X);

        At = T*sys.A/T;
        Bt = T*sys.B;
        Ct = sys.C/T;

        Q = eye(size(sys.A));
        E = Q;
        J = 0.5*(At - At');
        R = -0.5*(At + At');
        G = 0.5*(Ct' + Bt);
        P = 0.5*(Ct' - Bt);

    case 'implicit'
        E = X;
        Ai = E*sys.A;
        Bi = E*sys.B;

        Q = eye(size(X));
        J = 0.5*(Ai - Ai');
        R = -0.5*(Ai + Ai');
        G = 0.5*(sys.C' + Bi);
        P = 0.5*(sys.C' - Bi);

        % End of switch Opts.rep
end

% Construct S and N
S = 0.5*(D + D');
N = 0.5*(D - D');

% Correct P if D is singular
if Opts.Dsingular
    S = zeros(size(S));
    N = zeros(size(N));
    P = zeros(size(P));
end

% Construct pH system
sysPH = phsRed(J, R, Q, G, E, P, S, N, Opts.phs);
sysPH.method = Opts.method;
sysPH.parameters = Opts;

end

%% Auxilliary functions
function [sys,D,X,availableX,Opts] = parseInputs(sys,varargin)
% Check input
if nargin > 1 % More than sys is provided
    % sys and Opts provided
    if nargin==2 && isa(varargin{1},'struct')
        Opts = varargin{1};
        availableX = false;
        X = [];
        % sys and X provided
    elseif nargin==2
        X = varargin{1};
        Opts=struct();
        availableX = true;
        % sys, X and Opts provided
    else
        X = varargin{1};
        Opts = varargin{2};
        availableX = true;
    end
    % Only sys provided
else
    Opts = struct();
    availableX = false;
    X = [];
end

% Check type of sys
if isa(sys,'phs') || isa(sys,'phsRed')
    error('System is already in pH representation.')
elseif ~(isa(sys,'sss') || isa(sys,'ssRed') || isa(sys,'ss'))
    error('Input system is not of type sss, ssRed or ss.')
end

% Check if sys is DAE
if svds(sparse(sys.A),1,'smallest') < 1e-12
    error('ss2phs does not work with DAE systems.')
end

% Check if E exists
if ~isempty(sys.E)
    % Make sys is explicit
    sys.A = sys.E\sys.A;
    sys.B = sys.E\sys.B;
    sys.E = eye(size(sys.A));
end

% Set default values
OptsAdmissible.tol = 1e-12;
OptsAdmissible.representation = {'implicit','standard','scaled'};
OptsAdmissible.method = @ss2phs;
OptsAdmissible.prl = {'cvx','icare','mmess','yalmip'};
OptsAdmissible.samPassive.plot = {false,true};
OptsAdmissible.phs.inputTolerance = 1e-8;
OptsAdmissible.yalmip.solver = 'sedumi';

% MMESS options
OptsAdmissible.mmess.norm = 'fro';
OptsAdmissible.mmess.ri.riccati_solver = 'radi';
OptsAdmissible.mmess.ri.lqg_solver = 'radi';
OptsAdmissible.mmess.ri.maxiter = 500;
OptsAdmissible.mmess.ri.res_tol = 1e-15;
OptsAdmissible.mmess.ri.rel_diff_tol = 1e-15;
OptsAdmissible.mmess.ri.trunc_tol = 1e-15;
OptsAdmissible.mmess.ri.compres_tol = 1e-15;
OptsAdmissible.mmess.ri.info = 0;
OptsAdmissible.mmess.ri.Z0 = [];
OptsAdmissible.mmess.ri.store_lqg = 0;
OptsAdmissible.mmess.ri.store_solfac = 0;
OptsAdmissible.mmess.ri.trunc_info = 0;

OptsAdmissible.mmess.shifts.num_desired = 10;
OptsAdmissible.mmess.shifts.method = 'projection';

OptsAdmissible.mmess.adi.maxiter = 500;
OptsAdmissible.mmess.adi.res_tol = 1e-15;
OptsAdmissible.mmess.adi.rel_diff_tol = 1e-15;
OptsAdmissible.mmess.adi.info = 0;
OptsAdmissible.mmess.adi.compute_sol_fac = 0;
OptsAdmissible.mmess.adi.accumulateK = 0;
OptsAdmissible.mmess.adi.accumulateDeltaK = 0;

OptsAdmissible.mmess.radi.maxiter = 500;
OptsAdmissible.mmess.radi.res_tol = 1e-15;
OptsAdmissible.mmess.radi.rel_diff_tol = 1e-15;
OptsAdmissible.mmess.radi.info = 0;

OptsAdmissible.mmess.oper = 'default';

Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

% Check D+D' > 0
Opts.Dsingular = false;
ev = eig(sys.D+sys.D');
if all(ev > 0)
    % Positive Definite
    D = sys.D;
else
    D = sys.D + Opts.tol*eye(size(sys.D));
    Opts.Dsingular = true;
end

% Check third party software
switch Opts.prl
    case 'cvx'
        thirdPartyCheck('CVX');
    case 'yalmip'
        thirdPartyCheck('YALMIP');
    case 'mmess'
        thirdPartyCheck('M.E.S.S.');
end

end