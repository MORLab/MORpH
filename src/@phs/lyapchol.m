function [R, S] = lyapchol(sys)
% LYAPCHOL - Computes the Cholesky factors R and S of the pH system's 
%            associated Lyapunov equations
%
% Syntax:
%   R = lyapchol(sys)
%   [R, S] = lyapchol(sys)
%   R = sys.lyapchol
%   [R, S] = sys.lyapchol
%
% Description:
%       Matrix R: X = R'*R
%       For explicit systems, the considered Lyapunov equation is:
%           (J-R)*Q*X + X*Q'*(J-R)' + (G-P)*(G-P)' = 0
%
%       For implicit systems, the considered Lyapunov equation is:
%           (J-R)*Q*X*E' + E*X*Q'*(J-R)' + (G-P)*(G-P)' = 0
%
%       Matrix S (optional output): Y = S'*S
%       For explicit systems, the considered Lyapunov equation is:
%           Y*(J-R)*Q + Q'*(J-R)'*Y + Q'*(G+P)*(G+P)'*Q = 0
%
%       For implicit systems, the considered Lyapunov equation is:
%           E'*Y*(J-R)*Q + Q'*(J-R)'*Y*E + Q'*(G+P)*(G+P)'*Q = 0
%
%       Note that S is not computed if it is not requested --> time saving
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs-object
%
% Output Arguments:
%       - R:    Factorization of X = R'*R of reachability Gramian X
%       - S:    Factorization of Y = S'*S of observability Gramian Y
%
% See Also:
%       lyapchol, phs
%
% References:
%       [1] MATLAB, lyapchol
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

A = (sys.J - sys.R)*sys.Q;
B = (sys.G - sys.P);
C = (sys.G + sys.P)'*sys.Q;
E = sys.E;

% Reachability Gramian
if ~exist("mess_lyap.m", 'file')
    if ~sys.isImplicit
        % Explicit system
        R = lyapchol(A,B);
    else
        % Implicit system
        R = lyapchol(A,B,E);
    end
else % Installation of MESS Toolbox found
    if ~sys.isImplicit
        % Explicit system
        R = mess_lyap(A, B)';
    else
        % Implicit system
        R = mess_lyap(A, B, [], [], E)';
    end
end

if nargout > 1  % Observability Gramian
    if ~exist("mess_lyap.m", 'file')
        if ~sys.isImplicit
            % Explicit system
            S = lyapchol(A',C');
        else
            % Implicit system
            S = lyapchol(A',C',E');
        end
    else % Installation of MESS Toolbox found
        if ~sys.isImplicit
            % Explicit system
            S = mess_lyap(A',C')';
        else
            % Implicit system
            S = mess_lyap(A', C', [], [], E')';
        end
    end
end

end