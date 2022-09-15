function sysr = prbt(sys, varargin)
% PRBT - Obtaining a passive reduced order system by positive real
%        balancing
% Syntax:
%   sysr = PRBT(sys)
%   sysr = PRBT(sys, r)
%   sysr = PRBT(sys, Opts)
%   sysr = PRBT(sys, r, Opts)
%
% Description:
%       sysr = prbt(sys, r, Opts) returns a passive reduced system with
%       given order r by Riccati or mixed Gramian balancing [1].
%
%       If r equals the system order, the ssRed object sysr is a balanced
%       realization.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:          sss object, containing passive LTI system
%       *Optional Input Arguments:*
%       - r:            desired reduced order
%       - Opts:  structure with execution parameters
%           - .truncation               Define how the model will be truncated
%                                       [{'redOrder'} / 'truncTol' / 'userIn' / 'errBound']
%           - .redErr                   Maximum H-infinity error
%                                       [{1e-2} / positive double]
%           - .truncTol                 Truncation tolerance for modified
%                                       Hankel singular values relative to the greatest
%                                       [{1e-2} / positive double]
%           - .hsvTol                   Tolerance for modified Hankel singular values
%                                       [{1e-15} / positive double]
%           - .tol                      Set artificial D=tol*eye(m) if D=0
%                                       [{1e-12} / positive double]
%           - .useCholFactors         	Use Cholseky factorization instead of PRARE solution
%                                       [{false} / true]
%           - .checkPassivity:          Check passivity of sys
%                                       [{false} / true]
%           - .mixedGramian             Define which balancing is used
%                                       [{'standard'} / 'mixedControl' / 'mixedObserve']
%           - .lyap                   	Method to solve lyapunov equation
%                                       [{'auto'} / 'mmess' / 'lyapchol']
%           - .are                      Method to solve ARE equation
%                                       [{'auto'} / 'mmess' / 'icare']
%           - .makePH                   Define if sysr should be pH or not
%                                       [{false} / true]
%           - .areOpts.*                Options that will be passed to
%                                       mess_lrri.
%                                       Please refer to documentation of the respective
%                                       algorithm.
%           - .lyapOpts.*               Options that will be passed to
%                                       mess_lradi. Please refer to documentation of the
%                                       respective algorithm.
%           - .ss2phs.*                 Options that will be passed to
%                                       ss2phs. Please refer to documentation of the respective
%                                       algorithm.
%           - .phs.*                    Options that will be passed to phs.
%                                       Please refer to documentation of the respective
%                                       algorithm.
%           - .samPassive.*             Options that will be passed to samPassive.
%                                       Please refer to documentation of the respective
%                                       algorithm.
%
% Output Arguments:
%       - sysr:     ssRed (phsRed) object, containing reduced LTI (port-Hamiltonian) system
%
% See Also:
%       phs, icare, mess_lradi, mess_lrri, lyapchol, ss2phs, samPassive,
%       demo_prbt
%
%
% References:
%       [1] K. Unneland, P. Van Dooren, and O. Egeland. A Novel Scheme for Positive
%           Real Balanced Truncation. In 2007 American Control Conference, 
%           pages 947–952, 2007.
%       [2] J. Phillips, L. Daniel, and L. M. Silveira. Guaranteed passive balancing
%           transformations for model order reduction. In Proceedings 2002 
%           Design Automation Conference (IEEE Cat. No.02CH37324), pages 52–57, 2002.
%       [3] U. Desai and D. Pal. "A transformation approach to stochastic 
%           model reduction." In: IEEE Transactions on Automatic Control 29.12 (1984), 
%           pp. 1097–1100.
%       [4] J. Saak, M. Köhler, and P. Benner. M-M.E.S.S.-2.1 – The
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
%               If you use the software package 'M-M.E.S.S.' [4], please refer
%               to its specific license.
%-----------------------------------------------------------------------

%% Parse inputs
[sys,r,D,Opts] = parseInputs(sys,varargin{:});

%% Solve Riccati Equations
if strcmp(Opts.are, 'mmess')
    % Setup eqn struct for M.-M.E.S.S.
    R       = chol(D+D');
    Rinv    = R\eye(size(D));

    % eqn = struct('A_', sys.A, 'B1', sys.B*Rinv, 'C1', Rinv'*sys.C, 'U', -sys.B*Rinv,...
    %     'V', (Rinv'*sys.C)', 'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
    %     'haveUV', 1);
    eqn = struct('A_', sys.A-sys.B*Rinv*(Rinv'*sys.C), 'B1', sys.B*Rinv, 'C1', Rinv'*sys.C,...
        'B2', zeros(size(sys.B)), 'C2', zeros(size(sys.C)), 'haveE', 0,...
        'haveUV', 0);

    % Operator manager for M.-M.E.S.S.
    oper = operatormanager(Opts.areOpts.oper);

    if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedObserve')
        % Input energy/required supply gramian R [1]
        eqn.type = 'N';
        solR = mess_lrri(eqn, Opts.areOpts, oper);

        X = solR.Z*solR.Z';
        res = sys.A*X+X*sys.A'+(X*sys.C'-sys.B)/(D+D')*(X*sys.C'-sys.B)';
        err = max(abs(res(:)));

        if err > 1e-8
            warning('Computational issues with solution of algebraic Riccati equation.\n');
            warning('Error of algebraic Riccati equation is %.8f.\n',err);
        end
    end

    if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedControl')
        % Output energy/available storage Gramian O [1]
        eqn.type = 'T';
        solO = mess_lrri(eqn, Opts.areOpts, oper);

        X = solO.Z*solO.Z';
        res = sys.A'*X+X*sys.A+(X*sys.B-sys.C')/(D+D')*(X*sys.B-sys.C')';
        err = max(abs(res(:)));

        if err > 1e-8
            warning('Computational issues with solution of algebraic Riccati equation.\n');
            warning('Error of algebraic Riccati equation is %.8f.\n',err);
        end
    end

    if Opts.useCholFactors % Use Cholesky factors of solutions to positive real Riccati equations
        try
            % Required supply R
            if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedObserve')
                R = solR.Z*solR.Z';
                LR = chol(R, 'lower');
            end

            % Available storage O
            if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedControl')
                O = solO.Z*solO.Z';
                LO = chol(O, 'lower');
            end
        catch
            error('Gramian not positive definite. System is probably not minimal.')
        end

    else % Use solutions of positive real Riccati equations
        % Required supply R
        if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedObserve')
            LR = solR.Z;
        end

        % Available storage O
        if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedControl')
            LO = solO.Z;
        end
    end

else
    if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedObserve')
        % Input energy/required supply Gramian R [1]
        R = icare(full(sys.A'),full(sys.C'),0,full(-D-D'),full(-sys.B));

        [UR,DR] = eig(0.5*(R+R'));
        LR = (UR*sqrt(DR));

        res = sys.A*R+R*sys.A'+(R*sys.C'-sys.B)/(D+D')*(R*sys.C'-sys.B)';
        err = max(abs(res(:)));

        if err > 1e-8
            warning('Computational issues with solution of algebraic Riccati equation.\n');
            warning('Error of algebraic Riccati equation is %.8f.\n',err);
        end
    end
    if strcmp(Opts.mixedGramian, 'standard') || strcmp(Opts.mixedGramian, 'mixedControl')
        % Output energy/available storage Gramian O [1]
        O = icare(full(sys.A),full(sys.B),0,full(-D-D'),full(-sys.C'));

        [UO,DO] = eig(0.5*(O+O'));
        LO = (UO*sqrt(DO));

        res = sys.A'*O+O*sys.A+(O*sys.B-sys.C')/(D+D')*(O*sys.B-sys.C')';
        err = max(abs(res(:)));

        if err > 1e-8
            warning('Computational issues with solution of algebraic Riccati equation.\n');
            warning('Error of algebraic Riccati equation is %.8f.\n',err);
        end
    end
end

%% Solve Lyapunov equations
if ~strcmp(Opts.mixedGramian,'standard')
    if strcmp(Opts.lyap,'mmess')
        % Setup eqn struct for M.-M.E.S.S.
        lyapEqn=struct('A_', sys.A, 'B', sys.B, 'C', sys.C, 'haveE', 0);

        % Operator manager for M.-M.E.S.S.
        oper = operatormanager(Opts.areOpts.oper);

        % Get adi shifts
        Opts.lyapOpts.shifts.p = mess_para(lyapEqn, Opts.lyapOpts, oper);

        % Select which Lyapunov equation has to be solved
        if strcmp(Opts.mixedGramian, 'mixedControl')
            lyapEqn.type = 'N';
        elseif strcmp(Opts.mixedGramian, 'mixedObserve')
            lyapEqn.type = 'T';
        end

        % Low rank adi
        sol = mess_lradi(lyapEqn, Opts.lyapOpts, oper);
        LG = sol.Z;

        X = sol.Z*sol.Z';
    else
        if strcmp(Opts.mixedGramian, 'mixedControl')
            LG = lyapchol(sys.A,sys.B*sys.B')';
        elseif strcmp(Opts.mixedGramian, 'mixedObserve')
            LG = lyapchol(sys.A',sys.C'*sys.C)';
        end

        X = LG*LG';
    end

    % Check solution accuracy
    if strcmp(Opts.mixedGramian, 'mixedControl')
        res = sys.A*X + X*sys.A' + sys.B*sys.B';
    elseif strcmp(Opts.mixedGramian, 'mixedObserve')
        res = sys.A'*X + X*sys.A + sys.C'*sys.C;
    end

    err = max(abs(res(:)));

    if err > 1e-8
        warning('Computational issues with solution of Lyapunov equation.\n');
        warning('Error of equation is %.8f.\n',err);
    end
    % End of if ~strcmp(Opts.mixedGramian,'standard')
end

%% Singular value decomposition for balanced realization
if strcmp(Opts.mixedGramian, 'standard')
    [U,Sigma,V] = svd(LO'*LR, 0);
    hsv = diag(Sigma);

elseif strcmp(Opts.mixedGramian, 'mixedControl')
    [U,Sigma,V] = svd(LO'*LG, 0);
    hsv = diag(Sigma);

else % strcmp(Opts.mixedGramian, 'mixedObserve')
    [U,Sigma,V] = svd(LG'*LR, 0);
    hsv = diag(Sigma);
end

%% Determine reduction order
% Max reduced order
rmax=size(V, 1);
switch Opts.truncation
    case 'truncTol'
        r = find(diag(Sigma)/Sigma(1,1) >= Opts.truncTol, 1, 'last');

    case 'redOrder'
        % Check if r does exceed limits
        if r > rmax
            warning(['Reduced order exceeds maximum. r has been changed to ',...
                'rmax = ', num2str(rmax,'%d'), '.'])
            r = rmax;
        end

    case 'userIn'
        % Maximum reduced order
        rmax = min([sum(hsv >= Opts.hsvTol*hsv(1)), rmax]);

        % Setup figure and plot modified Hankel singular values
        h=figure;
        bar(1:rmax,abs(hsv(1:rmax)), 'r');
        %semilogy(1:rmax, hsv(1:rmax), 'ro')   % Alternative plot style

        % Axis properties
        title('Modified Hankel Singular Values');
        xlabel('Order');
        ylabel({'Relative hsv decay';sprintf('abs(hsv/hsv(1)) with hsv(1)=%.4d', hsv(1))});
        set(gca,'YScale','log');
        set(gca, 'YLim', [-Inf;1.5]);
        set(gca, 'XLim', [0; rmax]);

        % Get reduced order from user
        prompt=['Please enter the desired order (0<= r <=', num2str(rmax, '%d)'),': '];
        r = input(prompt);

        % Close figure
        if ishandle(h)
            close Figure 1;
        end

        % Check if r is not negative and an integer
        if r<0 || round(r)~=r
            error('Invalid reduction order.');
        end

        % Check if r does exceed limits
        if r > rmax
            warning(['Reduced order exceeds maximum. r has been changed to ',...
                'rmax = ', num2str(rmax,'%d'), '.'])
            r = rmax;
        end

    case 'errBound'
        % TODO: Reduction based on relative error bound [Benner (2004), Reis (2010)]
        error('Reduction via error bound is not possible yet.')

        % End of switch Opts.truncation
end

%% Truncation
U = U(:,1:r);
V = V(:,1:r);
Sigma = Sigma(1:r,1:r);

% Transformation matrices
if strcmp(Opts.mixedGramian, 'standard')
    Tinv = Sigma^-0.5*U'*LO';
    T = LR*V*Sigma^-0.5;

elseif strcmp(Opts.mixedGramian, 'mixedControl')
    Tinv = Sigma^-0.5*U'*LO';
    T = LG*V*Sigma^-0.5;

else % strcmp(Opts.mixedGramian, 'mixedObserve')
    Tinv = Sigma^-0.5*U'*LG';
    T = LR*V*Sigma^-0.5;
end

%% Create reduced order model
Ar = Tinv*sys.A*T;
Br = Tinv*sys.B;
Cr = sys.C*T;

sysr = ssRed(Ar,Br,Cr,sys.D,eye(size(Ar)),'prbt',Opts);

% Make reduced pH system
if Opts.makePH
    sysr = ss2phs(sysr, Opts.ss2phs);
    sysr.method = @prbt;
    sysr.parameters = Opts;
end

% End of function
end

%% Auxilliary functions
function [sys,r,D,Opts] = parseInputs(sys,varargin)
% Check input for r and Opts
if nargin>1
    % sys and Opts provided
    if nargin == 2 && ~isa(varargin{1}, 'double')
        Opts=varargin{1};

        % sys, r and possibly Opts provided
    else
        r = varargin{1};

        % Opts provided
        if nargin == 3
            Opts = varargin{2};
        else
            Opts = struct();
        end
    end
else
    Opts = struct();
end

% Check admissible option values
OptsAdmissible.useCholFactors = {false,true};   % Use Cholseky factorization instead of PRARE solution
OptsAdmissible.makePH = {false,true};           % Define if sysr should be pH or not
OptsAdmissible.ss2phs = struct();               % Make sure options struct for ss2phs exists
OptsAdmissible.redErr = 1e-2;                   % Error bound for H-infinity norm error
OptsAdmissible.truncTol = 1e-2;                 % Truncation tolerance (redOrder) for modified Hankel singular values
OptsAdmissible.hsvTol = 1e-15;                  % Tolerance for modified Hankel singular values
OptsAdmissible.phs = struct();                  % Make sure options struct for phs exists
OptsAdmissible.samPassive.plot = false;         % Disable plot in sampassive
OptsAdmissible.tol = 1e-12;                     % Set D=tol*eye(m) for D=0
OptsAdmissible.checkPassivity = {false,true};
OptsAdmissible.mixedGramian = {'standard','mixedControl','mixedObserve'};   % Define which balancing is used
OptsAdmissible.lyap = {'auto','mmess','lyapchol'};
OptsAdmissible.are = {'auto','mmess','icare'};
OptsAdmissible.truncation = {'redOrder','truncTol','userIn','errBound'};

% Options for MMESS ARE
OptsAdmissible.areOpts.ri.riccati_solver = 'radi';
OptsAdmissible.areOpts.ri.lqg_solver = 'radi';
OptsAdmissible.areOpts.ri.maxiter = 500;
OptsAdmissible.areOpts.ri.res_tol =  1e-15;
OptsAdmissible.areOpts.ri.rel_diff_tol = 1e-15;
OptsAdmissible.areOpts.ri.trunc_tol = 1e-15;
OptsAdmissible.areOpts.ri.compres_tol = 1e-15;
OptsAdmissible.areOpts.ri.info = 0;
OptsAdmissible.areOpts.ri.Z0 = [];
OptsAdmissible.areOpts.ri.store_lqg = 0;
OptsAdmissible.areOpts.ri.store_solfac = 0;
OptsAdmissible.areOpts.ri.trunc_info = 0;
OptsAdmissible.areOpts.norm = {'fro',2};
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

% Options for MMESS lyap
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

% Check for pH system
isPHS = 0;
if isa(sys, 'phs')
    sys = phs2sss(sys);
    isPHS = 1;
end

% Check for passive system
if ~isPHS && Opts.checkPassivity
    passive = samPassive(sys,Opts.samPassive);
else
    passive = true;
end

if ~passive
    error('Positive real balanced truncation is only applicable to passive systems.')
end

% Check for D+D' > 0 and correct if necessary
ev = eig(sys.D+sys.D');
if all(ev>0)
    % Positive Definite
    D = sys.D;
else
    D = sys.D + Opts.tol*eye(size(sys.D));
end

% Check if sys is DAE
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('prbt currently only supports ODE systems.')
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

% Check input according to Opts.truncation
switch Opts.truncation
    case 'redOrder'
        if ~exist('r', 'var')
            error('Reduced order r must be provided.')
        end
        if r <= 0 || round(r) ~= r
            error('Reduced order r must be a positive integer.')
        end
    case 'truncTol'
        if Opts.truncTol <= 0
            error('Truncation tolerance must be greater than zero.')
        end
    case 'errBound'
        error('Reduction via error bound is not possible yet.')
end

% Select lyap solver if 'auto'
if strcmp(Opts.lyap,'auto')
    if sys.n > 1000
        Opts.lyap = 'mmess';
    else
        Opts.lyap = 'lyapchol';
    end
end
% Select are solver if 'auto'
if strcmp(Opts.are,'auto')
    if sys.n > 1000
        Opts.are = 'mmess';
    else
        Opts.are = 'icare';
    end
end
% Check third party software
if strcmp(Opts.lyap,'mmess') || strcmp(Opts.are,'mmess')
    thirdPartyCheck('M.E.S.S.');
end

end
