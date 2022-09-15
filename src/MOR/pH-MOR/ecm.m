function [redSys, hsv] = ecm(sys, varargin)
% ECM - computes the reduced port-Hamiltonian system model redSys for the
%       full-order model sys using the effort-constraint method [1]
%
% Syntax:
%   redSys = ECM(sys)
%   redSys = ECM(sys, redOrder)
%   redSys = ECM(sys, redOrder, Opts)
%   redSys = ECM(sys, Opts)
%   [redSys, sigma] = ECM(sys, ...)
%
% Description:
%       redSys = ecm(sys, redOrder) computes the reduced
%       port-Hamiltonian system redSys.
%       The system will be transformed to a balanced co-energy
%       representation first. Subsequently, the negligible directions (in
%       energy coordinates) will be truncated. [1,2]
%
%       Note that the algorithm applies some sort of symmetry correction to
%       the reduced matrices J, R, and Q to prevent phs input validation
%       failures.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       - redOrder: desired order of the reduced system
%       *Optional Input Arguments:*
%       - Opts:     structure with execution parameters
%           - .symmetryCorrection:
%                       Decide wheter matrices J, R, Q will be 'corrected'
%                       in terms of symmetry.
%                       [{true} / false]
%           - .tolSingularValues: Value below which the algorithm will
%                       advise to neglect Hankel singular values.
%                       [{1e-12} / positive double]
%           - .phs.*:   Options that will be passed on to the 'phs' class.
%
% Output Arguments:
%       - redSys: 	reduced order port-Hamiltonian model (phsRed object)
%       - hsv:      Hankel singular values of the original system
%
% See Also:
%       demo_ecm, phs
%
% References:
%       [1] Polyuga (2010), "Model Reduction of Port-Hamiltonian Systems",
%           Dissertation, 2010.
%       [2] R. V. Polyuga and A. van der Schaft. "Effort- and Flow-Constraint 
%           Reduction Methods for Structure Preserving Model Reduction of 
%           Port-Hamiltonian Systems.” In: Systems Control Lett. 61.3 (2012), 
%           pp. 412–421.
%       [3] A. C. Antoulas. Approximation of Large-Scale Dynamical Systems, 
%           Adv. Des. Control 6. SIAM, Philadelphia, 2005.
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Julius Durmann, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Parse inputs
narginchk(1,3);

if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end
OptsAdmissible.phs = struct();
OptsAdmissible.tolSingularValues = 1e-12;
OptsAdmissible.symmetryCorrection = true;
Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

if ~isempty(varargin)
    redOrder = varargin{1};
else
    redOrder = -1;
end

% Catch incompatible DAE systems
if sys.isDAE
    error('ecm currently only supports pHODE systems.')
end
% Make implicit systems explicit
if sys.isImplicit
    warning('Transforming descriptor system (implicit) to explicit system... Computation of matrix inverse is performed!')
    warning('off', 'phs:phs:ChangedProperty')
    sys = sys.makeExplicit();
    warning('on', 'phs:phs:ChangedProperty')
end

%% Balancing (see [3])
% Calculate state-transformation T: e_bal = T*e for the co-energy
% system. This transformation balances the system for the co-energy
% coordinates.

try
    A = (sys.J - sys.R)*sys.Q;
    B = (sys.G - sys.P);
    C = (sys.G + sys.P)'*sys.Q;
    S = lyapchol(A,B);         % returns W_c = S'*S, W_c: controllabiltiy Gramian
    R = lyapchol(A',C');       % returns W_o = R'*R, W_o: observability Gramian
    S = S*sys.Q';  % Lyapunov equations need to be solved for the co-energy coordinate system.
    R = R/sys.Q;   % The standard solution is for the energy coordinate system.
catch exception
    % Lyapchol errored
    msg = ['Computation of the Gramians (using lyapchol) errored. '...
        'The realization of the system is probably not minimal and '...
        'the Gramians are not positive definite.'];
    causeException = MException('MORpH:ecm:minimalRealizationError',msg);
    exception = addCause(exception,causeException);
    rethrow(exception)
end

[U, SIGMA, V] = svd(S*R');
hsv = diag(SIGMA);    % Vector of singular values (for element-wise computation)
nValid = find(hsv > Opts.tolSingularValues, 1, 'last');

T = (hsv.^(-0.5)).*V'*R;      % state transformation matrix T
T_inv = S'*U.*(hsv.^(-0.5))';  % state transformation matrix T^(-1)

J_bal = T_inv'*sys.J*T_inv;
R_bal = T_inv'*sys.R*T_inv;
Q_bal = T*sys.Q*T';
G_bal = T_inv'*sys.G;
P_bal = T_inv'*sys.P;

%% Request user input of reduced order
% Ask user if no reduced order has been passed
if redOrder < 0
    f = figure();
    semilogy(1:length(hsv), hsv, '*')
    title('Hankel singular values')
    redOrder = str2double(input(['Please enter a reduced order based on the Hankel\n'...
        'singular values of the system. '], 's'));
    if ishandle(f), close(f); end
end
% Ask user (again) if minimal order is less than reduced order
if nValid < redOrder
    f = figure();
    semilogy(1:length(hsv), hsv, '*')
    title('Hankel singular values')
    redOrder = str2double(input(['The reduced order you requested (',...
        num2str(redOrder), ') \nis greater than the minimal realization '...
        'order (', num2str(nValid), ') \n'...
        'of the system.\n'...
        'With which order would you like to proceed? '], 's'));
    if ishandle(f), close(f); end
end

%% Reduction

% Truncation according to [1,2]
J_red = J_bal(1:redOrder,1:redOrder);
R_red = R_bal(1:redOrder,1:redOrder);
Q_help = Q_bal(redOrder+1:end,redOrder+1:end) \ Q_bal(redOrder+1:end,1:redOrder);
Q_red = Q_bal(1:redOrder,1:redOrder) - Q_bal(1:redOrder, redOrder+1:end) * Q_help;
G_red = G_bal(1:redOrder, :);
P_red = P_bal(1:redOrder, :);

% Create reduced system
if Opts.symmetryCorrection
    [J_red,R_red,Q_red] = correctSymmetry(J_red,R_red,Q_red);
end

redSys = phsRed(J_red, R_red, Q_red, G_red, eye(size(Q_red)), P_red, sys.S, sys.N, Opts.phs);
redSys.parameters = Opts;
redSys.method = @ecm;

end

function [J_out, R_out, Q_out] = correctSymmetry(J, R, Q)
% Correct slightly unsymmetric matrices J, R so that
%   J_out' = -J_out,
%   R_out' = R_out
%   J_out - R_out = J - R
% The reason for this correction is that the phs input validation often
% fails because of missing symmetry properties. This function therefore
% makes ecm more robust.

J_out = 0.5*((J-R)-(J-R)');
R_out = -0.5*((J-R)+(J-R)');
if nargin == 3
    Q_out = 0.5*(Q+Q');
end

errJ = norm(J_out-J)/norm(J);
errR = norm(R_out-R)/norm(R);

msg = ['To avoid failure of phs.inputValidation, symmetry correction was '...
    'applied to the matrices J, R, and Q. This makes ecm '...
    'more robust. To turn off this feature, use Opts.symmetryCorrection = false.\n'];
if nargin == 3
    errQ = norm(Q_out-Q)/norm(Q);
    warningMsg = strcat(msg, 'By symmetry correction the following relative changes were made: \nerrJ = ', ...
        num2str(errJ),', errR = ', num2str(errR), ', errQ = ', num2str(errQ));
else
    warningMsg = strcat(msg, 'By symmetry correction the following relative changes were made: \nerrJ = ', ...
        num2str(errJ),', errR = ', num2str(errR));
end
warning('MORpH:ecm:correctSymmetry', warningMsg);

end
