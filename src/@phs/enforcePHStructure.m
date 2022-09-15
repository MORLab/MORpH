function [sys, changes] = enforcePHStructure(sys, proceedAnyway)
% ENFORCEPHSTRUCTURE - Modifies system matrices of sys such that the system
% fulfills the pH structural constraints
%
% Syntax:
%             sys = sys.ENFORCEPHSTRUCTURE
%   [sys,changes] = sys.ENFORCEPHSTRUCTURE(proceedAnyway)
%
% Description:
%       This function modifies the system matrices of sys (if
%       necessary) to make it a port-Hamiltonian system. Therefore, (skew-)
%       symmetry and positive (semi-) definiteness of the matrices is
%       checked and enforced.
%
%       Note: This function is not supposed to transform arbitrary systems
%       to PH models but is instead intended to correct minor deviations
%       from the desired structure which may occur, for instance, due to
%       model order reduction. If deviations are large, the returned model
%       will differ significantly from the original. The function is only
%       applicable to systems of small or medium size.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:	system with small violations of the pH structural
%               constraints
%       *Optional Input Arguments:*
%       - proceedAnyways:	Can be used to avoid user input prompt
%                           [{false} / true]
% Output Arguments:
%       - sys:      Corrected system which fulfills the pH structural
%                   constraints
%       - changes:  Summary of the changes being made
%
% Examples:
%       The code below creates a PH Mass-Spring-Damper model, adds some
%       deviations from the PH structure to it and runs the
%       enforcePHStructure algorithm.
%
%       sys = setup_MassSpringDamperSystem(50, 2, 1, 1);
%       sys.J(end,1)= sys.J(end,1) + 1;
%       sys.Q = sys.Q - 1e-2*eye(sys.dim);
%       [sys_corrected, changes] = sys.enforcePHStructure
%
% See Also:
%       phs
%
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

%% Check system size
narginchk(1,2)
% Catch large systems to avoid unintended large computations
if nargin < 2
    proceedAnyway = false;
end
warning('off','phs:phs:ChangedProperty')
if sys.dim > 5e1
    warning('MORpH:phs:enforcePHStructure:largeSystem',...
        'System dimension is large. Computation may take some time.');
    if ~proceedAnyway
        answer = input('Do you want to proceed? (y/n)    ','s');
        if ~ (strcmp(answer, 'Y') || strcmp(answer, 'y') || isempty(answer))
            fprintf(2,'Aborted computation ');
            return
        else
            disp('Continuing...');
        end
    end
end

sys_old = sys;

%% Make Q'*E symmetric positive semi-definite
QtE = sys.Q'*sys.E;
QtE = makePSD(full(QtE),1e-12);

if rank(full(sys.E)) == sys.dim
    sys.Q = (QtE/sys.E)';
elseif rank(full(sys.Q)) == sys.dim
    sys.E = (sys.Q')\QtE;
else
    warning("phs:enforcePHstructure:changedAboveTolerance",...
        "Matrix product Q'*E can only be corrected if either Q or E is invertible.")
end

%% Make Q'*J*Q skew-symmetric and Q'*R*Q symmetric
if ~issymmetric(sys.Q'*sys.J*sys.Q,'skew') || ~issymmetric(sys.Q'*sys.R*sys.Q)
    sys.J = 0.5*((sys.J-sys.R)-(sys.J-sys.R)'); % skew-symmetric
    sys.R = -0.5*((sys.J-sys.R)+(sys.J-sys.R)');% symmetric
end

%% Make N skew-symmetric and S symmetric
if ~issymmetric(sys.S) || ~issymmetric(sys.N,'skew')
    sys.S = 0.5*((sys.N+sys.S)+(sys.N+sys.S)'); % symmetric
    sys.N = 0.5*((sys.N+sys.S)-(sys.N+sys.S)'); % skew-symmetric
end

%% Make W = [Q'*R*Q, Q'*P; P'*Q, S] positive semi-definite
W = full([sys.Q'*sys.R*sys.Q, sys.Q'*sys.P; sys.P'*sys.Q, sys.S]);
W = makePSD(W,1e-12);
sys.R = sys.Q'\W(1:sys.dim,1:sys.dim)/sys.Q;
sys.P = sys.Q'\W(1:sys.dim,sys.dim+1:end);
sys.S = W(sys.dim+1:end,sys.dim+1:end);

%% Compute absolute and relative changes
matrices = {'J','R','Q','G','E','P','S','N'};
for i = 1:length(matrices)
    mat = matrices{i};
    changes.abs.(mat) = norm(full(sys.(mat)-sys_old.(mat)));
    changes.rel.(mat) = changes.abs.(mat)/norm(full(sys.(mat)));
end

%% Console output
if sys.Opts.verbose
    fprintf(" Absolute changes\t\t|\tRelative Changes\n");
    fprintf("-------------------------------------------\n");
end
for i = 1:length(matrices)
    mat = matrices{i};
    if sys.Opts.verbose
        fprintf(" %s-%s_old = %.3e\t|\t",mat,mat,changes.abs.(mat));
        fprintf("%s-%s_old/%s_old = %.3e\n",mat,mat,mat,changes.rel.(mat));
    end
    % Report large changes / Issue warning
    if  changes.rel.(mat) > sys.Opts.inputTolerance
        msg = strcat("Changed property ", mat, " above tolerance (", string(sys.Opts.inputTolerance), ")");
        warning("phs:enforcePHstructure:changedAboveTolerance", msg);
    end
end

warning('on','phs:phs:ChangedProperty')
try
    phs.inputValidation(sys);
catch
    warning("phs:enforcePHstructure:correctedSystemNotPH", ...
        "PH structure could not be enforced.");
end

end