function sys = mtimes(sys1,sys2)
% MTIMES - returns a system which is obtained by multiplying sys1 and sys2.
%          Can be called by sys1*sys2.
%
% Description:
%   u ------> [ sys2 ] --> [ sys1 ] ---------> y
%   This method is not structure preserving in general!
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1:	first system (phs object)
%       - sys2: second system (phs object)
%
% Output Arguments:
%       - sys:  sys1*sys2 (system chain) (ss object)
%
% See Also:
%      	phs
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

if (isnumeric(sys1) && isreal(sys1) && isscalar(sys1)) && isa(sys2,'phs')
    if sys1 >= 0
        [J, R, Q, G, E, P, S, N] = getMatrices(sys2);
        G = sqrt(sys1)*G;
        P = sqrt(sys1)*P;
        S = sys1*S;
        N = sys1*N;
        sys = phs(J,R,Q,G,E,P,S,N);
    else
        % Convert to ss
        sys1 = ss(sys1);
        sys2 = ss(sys2);
        sys = sys1*sys2;
    end

elseif isa(sys1,'phs') && (isnumeric(sys2) && isreal(sys2) && isscalar(sys2))
    if sys2 >= 0
        [J, R, Q, G, E, P, S, N] = getMatrices(sys1);
        G = sqrt(sys2)*G;
        P = sqrt(sys2)*P;
        S = sys2*S;
        N = sys2*N;
        sys = phs(J,R,Q,G,E,P,S,N);
    else
        % Convert to ss
        sys1 = ss(sys1);
        sys2 = ss(sys2);
        sys = sys1*sys2;
    end


elseif isa(sys1,'phs') && (isa(sys2,'phs') || isa(sys2,'ss') || isa(sys2,'sss'))
    % Convert to ss
    sys1 = ss(sys1);
    sys2 = ss(sys2);
    sys = sys1*sys2;

else
    error('phs:mtimes:wrongInput',...
        'phs objects may be multiplied by positive real scalars and phs, ss and sss objects.');
end