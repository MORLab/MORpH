function sysr = balPH(sys, varargin)
% balPH - obtains a reduced order port-Hamiltonian system by
%         balancing w.r.t. modified Gramians.
%
% Syntax:
%   sysr = balPH(sys)
%   sysr = balPH(sys, r)
%   sysr = balPH(sys, Opts)
%   sysr = balPH(sys, r, Opts)
%
% Description:
%       sysr = balPH(sys, r, Opts) returns a reduced PH system with
%       given order r by balanced truncation with a modified controllability
%       Lyapunov equation [1].
%
%       If r equals the system order, the phs-object sysr is the balanced
%       realization.
%
%       Used options, modified Hankel singular values and projection matrix
%       Wr are stored in the phs object sysr.
%
%       If a reduced order is passed to the function the reduced system
%       will be this size. In this case Opts.truncTol is ignored. If no reduced
%       order is provided the system is truncated according to Opts.truncTol
%       if the value is greater than zero.  The function will
%       truncate the system according to the modified Hankel singular values
%       above the truncation tolerance. The option can be avoided
%       by setting Opts.truncTol to zero. If done so, the modified
%       Hankel singular values will be plotted for the user to enter a desired
%       reduced order.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:          phs object, containing LTI port-Hamiltonian system
%       *Optional Input Arguments:*
%       - r:            desired reduced order
%       - Opts:  structure with execution parameters
%           - .truncTol                 truncation tolerance for modified Hankel singular values
%                                       [{0} / positive integer]
%           - .hsvTol                   tolerance for modified Hankel singular values
%                                       [{1e-15} / positive integer]
%           - .lyap                     solver for Lyapunov equation
%                                       [{'mmess'} / 'lyapchol']
%           - .mess_lradi.*             Other options that will be passed on to the
%                                       M-M.E.S.S. method mess_lradi
%                                       Please refer to documentation of the respective
%                                       algorithm.
%
% Output Arguments:
%       - sysr:         phsRed object, containing reduced LTI port-Hamiltonian system
%
% See Also:
%       phs, mess_lradi, ecm, demo_balPH
%
%
% References:
%       [1] T. Breiten, R. Morandin, and P. Schulze. Error bounds for port-
%           Hamiltonian model and controller reduction based on system balancing. 
%           Comput. Math. Appl., 2021.
%       [2] J. Saak, M. Köhler, and P. Benner. M-M.E.S.S.-2.1 – The
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
%               If you use the software package 'M-M.E.S.S.' [2], please refer
%               to its license.
%-----------------------------------------------------------------------

%% Parse inputs
[sys,r,Opts] = parseInputs(sys,varargin{:});

%% Solve Lyapunov equations and get Cholesky factors
if strcmp(Opts.lyap,'mmess')
    % Setup eqn struct for M.-M.E.S.S.
    eqn=struct('A_', (sys.J-sys.R)*sys.Q, 'B', (sys.G-sys.P), 'C', (sys.G+sys.P)'*sys.Q, 'type', 'T', 'haveE', 0);

    % Get operatormanager and norm for M.-M.E.S.S.
    oper = operatormanager(Opts.mess_lradi.oper);

    % Get adi shifts for M.-M.E.S.S.
    Opts.mess_lradi.shifts.p = mess_para(eqn, Opts.mess_lradi, oper);

    % Low rank adi M.-M.E.S.S
    sol = mess_lradi(eqn, Opts.mess_lradi, oper);
    LY = sol.Z;
else
    A = (sys.J-sys.R)*sys.Q;
    C = (sys.G+sys.P)'*sys.Q;
    LY = lyapchol(A',C');
    LY=LY';
end

% Compute modified controllability Gramian X - inverse of Q
X = sys.Q\speye(size(sys.Q));

% Compute Cholesky factorization of Q
LX = chol(X);

%% Singular value decomposition for balanced realization
[~, Sigma, V] = svd(LX*LY, 0);

% Modified Hankel singular values
hsv = diag(Sigma);

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
        % End of switch Opts.truncation
end

%% Truncation
V = V(:,1:r);
Sigma = Sigma(1:r,1:r);

% Set modified hsv up to tolerance - prohibit bad conditioning of Q
w = find(diag(Sigma) < Opts.hsvTol*hsv(1));
idx = sub2ind(size(Sigma), w, w);
Sigma(idx) = Opts.hsvTol*hsv(1);

%% Create reduced order system (see [1] for details)
% Compute projection matrix Wr
Wr = LY*V*Sigma^-0.5;

% Compute system matrices of ROM
Jr = Wr'*sys.J*Wr;
Rr = Wr'*sys.R*Wr;
Qr = Sigma\eye(size(Sigma));
Gr = Wr'*sys.G;
Pr = Wr'*sys.P;
Er = speye(size(Qr));

% Get ROM
sysr = phsRed(Jr, Rr, Qr, Gr, Er, Pr, sys.S, sys.N, Opts.phs);

% Set properties of sysr
sysr.parameters = Opts;
sysr.parameters.Wr = Wr;
sysr.parameters.hsv = hsv;
if strcmp(Opts.lyap,'mmess')
    sysr.parameters.lyap = sol;
end
sysr.method = @balPH;

end

%% Auxilliary functions
function [sys,r,Opts] = parseInputs(sys,varargin)
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
OptsAdmissible.truncTol = 0;        % Truncation tolerance (redOrder) for modified Hankel singular values
OptsAdmissible.hsvTol = 1e-15;      % Tolerance for modified Hankel singular values
OptsAdmissible.phs = struct();      % Make sure options struct for phs exists
OptsAdmissible.lyap = {'mmess','lyapchol'};

% mmess Options
Opts.mess_lradi.adi.maxiter = 500;
Opts.mess_lradi.adi.res_tol = 0;
Opts.mess_lradi.adi.rel_diff_tol = 0;
Opts.mess_lradi.adi.info = 0;
Opts.mess_lradi.shifts.num_desired = 25;
Opts.mess_lradi.shifts.num_Ritz = 50;
Opts.mess_lradi.shifts.num_hRitz = 25;
Opts.mess_lradi.shifts.b0 = ones(sys.dim,1);
Opts.mess_lradi.shifts.method = 'projection';
Opts.mess_lradi.shifts.info = 0;
Opts.mess_lradi.oper = 'default';
Opts.mess_lradi.norm = 'fro';

Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

% Check for DAE system
if sys.isDAE
    error('balPH currently only supports pHODE systems.');
end

% Make implicit systems explicit
if sys.isImplicit
    warning('Transforming descriptor system (implicit) to explicit system... Computation of matrix inverse is performed!')
    warning('off', 'phs:phs:ChangedProperty')
    sys = sys.makeExplicit();
    warning('on', 'phs:phs:ChangedProperty')
end

% Determine truncation method
if exist('r', 'var') || Opts.truncTol > 0
    if exist('r', 'var')
        Opts.truncation = 'redOrder';

        % Check if r is not negative and an integer
        if r < 0 || round(r) ~= r
            error('Invalid reduction order.');
        end

    else
        Opts.truncation = 'truncTol';
    end
else
    Opts.truncation = 'userIn';
end

% Check third party software M.E.S.S.
thirdPartyCheck('M.E.S.S.')

% Check for r
if ~exist('r', 'var')
    r = [];
end

end