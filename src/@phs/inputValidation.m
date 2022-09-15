function isPH = inputValidation(sys)
% INPUTVALIDATION - checks if the system sys fulfills the pH structural constraints
%
% Syntax:
%   isPH = phs.INPUTVALIDATION(sys)
%
% Description:
%       phs.inputValidation(sys) checks the properties of the matrices of
%       sys according to the necessary PH properties introduced in [1],
%       definition 5
%
% Input arguments:
%   - sys:                      phs object
% Output arguments:
%   - isPH:                     true: system is pH
%                               false: system is not pH (should
%                               normally throw error)
%
% See Also:
%       phs, isPositiveDefinite
%
% References:
%       [1] C. Beattie, V. Mehrmann, H. Xu, and H. Zwart. Port-Hamiltonian descriptor systems. Math.
%           Control Signals Systems, 30(17):1–27, 2018.
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

%% Check if dimensions are correct

sE = size(sys.E);
% E not square
if (sE(1) ~= sE(2))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of E.');
end
sJ = size(sys.J);
% J not square or dimension different to E
if (sJ(1) ~= sJ(2) || ~isequal(sJ, sE))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of J.');
end
sR = size(sys.R);
% R not square or dimension different to E
if ( (sR(1) ~= sR(2)) || ~isequal(sE, sR) )
    error('MORpH:phs:wrongInput', 'Please check the dimensions of R.');
end
sQ = size(sys.Q);
% Q not square or dimension different to E
if ( (sQ(1) ~= sQ(2)) || ~isequal(sE, sQ) )
    error('MORpH:phs:wrongInput', 'Please check the dimensions of Q.');
end
sG = size(sys.G);
% G has wrong state dimension compared to E
if (sG(1) ~= sE(1))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of B.');
end
% Different dimensions of P and G
if (~isequal(sG, size(sys.P)))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of P.');
end
% S has wrong dimensions
if (size(sys.S,1) ~= sG(2) || size(sys.S,2) ~= sG(2))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of S.');
end
% N has wrong dimensions
if (size(sys.N,1) ~= sG(2) || size(sys.N,2) ~= sG(2))
    error('MORpH:phs:wrongInput', 'Please check the dimensions of N.');
end

%% Check (skew-)symmetry conditions
% See [1], Definition 5
errorMessage = "\nConsider changing Opts.inputTolerance to avoid this error when caused by numerical inaccuracy.";

QtE = sys.Q'*sys.E;
relErrQtE = abs(QtE - QtE')/max(abs(QtE(:)));
if max(relErrQtE(:)) > sys.Opts.inputTolerance ...
        && max(abs(QtE(:))) > sys.Opts.inputTolerance
    errorMessage = strcat("Matrix product Q'*E is not symmetric:\n Maximum relative error = ", num2str(max(relErrQtE(:))),errorMessage);
    error('MORpH:phs:wrongInput', errorMessage);
end

QtJQ = sys.Q'*sys.J*sys.Q;
relErrQtJQ = abs(QtJQ + QtJQ')/max(abs(QtJQ(:)));
if max(relErrQtJQ(:)) > sys.Opts.inputTolerance ...
        && max(abs(QtJQ(:))) > sys.Opts.inputTolerance
    errorMessage = strcat("Matrix product Q'*J*Q is not skew-symmetric:\n Maximum relative error = ", num2str(max(relErrQtJQ(:))),errorMessage);
    error('MORpH:phs:wrongInput', errorMessage);
end

relErrN = abs(sys.N + sys.N')/max(abs(sys.N(:)));
if max(relErrN(:)) > sys.Opts.inputTolerance ...
        && max(abs(sys.N(:))) > sys.Opts.inputTolerance
    errorMessage = strcat("N is not skew-symmetric:\n Maximum relative error = ", num2str(max(relErrN(:))),errorMessage);
    error('MORpH:phs:wrongInput', errorMessage);
end

W = full([sys.Q'*sys.R*sys.Q, sys.Q'*sys.P; sys.P'*sys.Q, sys.S]);
relErrW = abs(W - W')/max(abs(W(:)));
if max(relErrW(:)) > sys.Opts.inputTolerance ...
        && max(abs(W(:))) > sys.Opts.inputTolerance
    errorMessage = strcat("Dissipation matrix W = [Q'*R*Q, Q'*P; P'*Q, S] is not symmetric:\n Maximum relative error = ", num2str(max(relErrW(:))),errorMessage);
    error('MORpH:phs:wrongInput', errorMessage);
end

%% Check positive semidefiniteness of Q'E and W
% See [1], Definition 5

if (~phs.isPositiveDefinite(sys.Q'*sys.E,sys.Opts.inputTolerance,1,sys.Opts))
    errorMessage = strcat("Matrix product Q'*E is not positive semidefinite (Hamiltonian function must be bounded from below).",errorMessage);
    error('MORpH:phs:wrongInput', errorMessage)
end

if (~phs.isPositiveDefinite(W,sys.Opts.inputTolerance, 1, sys.Opts))
    errorMessage = strcat("Dissipation matrix W = [Q'*R*Q, Q'*P; P'*Q, S] is not positive semidefinite",errorMessage);
    error('MORpH:phs:wrongInput', errorMessage)
end

isPH = true;

end % of inputValidation()