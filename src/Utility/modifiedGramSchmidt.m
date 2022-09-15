function V_orth = modifiedGramSchmidt(V)
% MODIFIEDGRAMSCHMIDT - Returns an orthonormal basis V_orth of the column
%                       vectors in V
%
% Syntax:
%   V_orth = modifiedGramSchmidt(V)
%
% Description:
%       This function applies the modified Gram-Schmidt algorithm as
%       described in [1].
%
% Input Arguments:
%       *Required Input Arguments:*
%       - V: Matrix of column vectors
%
% Output Arguments:
%       - V_orth: Matrix with orthonormal columns which spans the same
%                 space as V
%
% References:
%       [1] https://www.math.uci.edu/~ttrogdon/105A/html/Lecture23.html
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Julius Durmann
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

if size(V, 2) > size(V, 1)
    error('MORpH:modifiedGramSchmidt:DimensionsNotSuitable',...
        'Cannot orthogonalize column matrix with more columns than rows!')
end

V_orth = V;
tol = 1e-12; % Tolerance for detecting zero vectors (linear dependence)

for j = 1:size(V,2)
    % Normalize current vector
    n = norm(V_orth(:,j));
    if all(abs(V_orth(:,j)) < tol)  % Vector is close to zero
        warning('MORpH:modifiedGramSchmidt:LinearlyDependentColumn',...
            'Found linearly dependend column. Replacing it by zero vector.');
        V_orth(:,j) = zeros(size(V,1),1);
        continue
    else
        V_orth(:,j) = V_orth(:,j) / n;
    end
    % Orthogonalize all subsequent vectors with respect to the current
    % one
    for k = j+1:size(V,2)
        V_orth(:,k) = V_orth(:,k) - (V_orth(:,k)'*V_orth(:,j))*V_orth(:,j);
    end
end

end