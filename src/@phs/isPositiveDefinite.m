function result = isPositiveDefinite(A,tol,semi,varargin)
% ISPOSITIVEDEFINITE - Determines whether the symmetric matrix A is positive
%       (semi-) definite. Uses eig or eigs (depending on dimension of
%       matrix A) to determine the smallest real part of the eigenvalues of
%       the eigenvalues of A. If this value is smaller than zero - tol, the
%       function will return false, otherwise it will return true
%
% Syntax:
%   result = ISPOSITIVEDEFINITE(A, tol, 0)
%   result = ISPOSITIVEDEFINITE(A, tol, 1)
%
% Description:
%       result = isPositiveDefiniteHermitian(A,0,1) returns logical value
%       'result' which is 1 for positive definite A and 0 otherwise
%
%       result = isPositiveDefiniteHermitian(A, 'semi') returns
%       value 'result' which is 1 for positive semi-definite A and 0
%       otherwise
%
% Input Arguments:
%       *Required Input Arguments:*
%       - A:		Hermitian matrix
%       - tol:      Tolerance
%       - semi:     0: check for positive definiteness
%                   1: check for positive semi-definiteness
%       *Optional Input Arguments:*
%       - Opts:     For activating/deactivating verbose mode
%                   (default: Opts.verbose = true)
%
% Output Arguments:
%       - result:   logical value
%
%
% See also:
%       eig, eigs
%
% References:
%       [1] Documentation of eigs
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

%% Input parsing
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
else
    Opts.verbose = true;
end

if length(A) > 5e3
    use = 'eigs'; %'eig' or 'eigs'
else
    use = 'eig';
end

narginchk(3,4);
result = true;


%% Compute smallest eigenvalues
switch use
    case 'eig'
        x = min(real(eig(full(A))));

    case 'eigs'
        EigsOpts.maxit = 5e2;
        EigsOpts.tol = tol;
        % EigOpts.disp = 1;
        x = eigs(A,1,'smallestreal',EigsOpts);
end

%% Check conditions
if ~semi
    % Positive definite?
    if isnan(x) && Opts.verbose
        warning('phs:isPositiveDefiniteHermitian:eigsThrewNan',...
            ['Could not ensure that the matrix is positive definite (', use, ' did not converge).'])
    else
        if x <= 0 - tol
            result = false;
        elseif x <= 0 && Opts.verbose
            warning('phs:isPositiveDefiniteHermitian:definitenessUnsure',...
                ['Matrix might not be completely positive definite '...
                '(computed eigenvalue is less than or equal to zero), '...
                'but within tolerance.\n'...
                'Smallest computed eigenvalue: ', num2str(x)]);
        end
    end
else
    % Positive semidefinite?
    if isnan(x) && Opts.verbose
        warning('phs:isPositiveDefiniteHermitian:eigsThrewNan',...
            ['Could not ensure that the matrix is positive semidefinite (', use, ' did not converge).'])
    else
        if x < 0 - tol
            result = false;
        elseif x < 0 && Opts.verbose
            warning('phs:isPositiveDefiniteHermitian:semiDefinitenessUnsure',...
                ['Matrix might not be completely positive semidefinite '...
                '(computed eigenvalue is less than zero), but within '...
                'tolerance.\n'...
                'Smallest computed eigenvalue: ', num2str(x)]);
        end
    end
end

if result == false && Opts.verbose
    warning(['Smallest computed eigenvalue: ', num2str(x)]);
end

end