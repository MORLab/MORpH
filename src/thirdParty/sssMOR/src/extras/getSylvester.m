function [R, B_, S] = getSylvester(sys,sysr,V,type)
% GETSYLVESTER - Get matrices of Sylvester's equation for Krylov subspaces
%
% Syntax: 
%       [Rv, B_, Sv] = GETSYLVESTER(sys, sysr, V)
%       [Lw, C_, Sw] = GETSYLVESTER(sys, sysr, W, 'W')
%
% Description:
%       Given the full order model sys and a reduced order model sysr
%       obtained through a projection with V being an input Krylov
%       subspace, this function reconstructs the matrices of the Sylvester
%       equation for the Krylov subspace
%
%           $A V - E V S_v - B R_v = 0 \quad          (1)$
%
%           $B\_ = (I - E V(W^T E V)^{-1}W^T) B \quad          (2)$
%
%           $A V - E V E_r^{-1} A_r - B\_ R_v = 0\quad   (3)$
%
%       If the type is set to 'W', then the matrices are given for the
%       output Krylov Sylvester equation
%
%           $A^T W - E^T W S_w^T - C^T L_w = 0 \quad   (4)$
%
%           $C\_ = C (I - V(W^T E V)^{-1}W^T E) \quad          (5)$
%
%           $A^T W - E^T V E_r^{-T} A_r^T - C\_^T L_w = 0 \quad  (6)$
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:      full order model
%       -sysr:     reduced order model
%       -V:        input Krylov subspace
%       *Optional Input Arguments:*
%       -type:     specifies if V spans an input (def) or output Krylov subspace 
%                  [{'V'} / 'W'] 
%
% Output Arguments: 
%       -R,S:      matrices of Sylvester equation (1) [Rv,Sv] or (4) [Lw,Sw]
%       -B_:       matrix of Sylvester equation (2) [B_] or (5) [C_]
% 
% Examples:
%       This code computes the input Krylov subspace for a benchmark model
%       and reconstructs the matrices of the corresponding Sylvester
%       equation
%
%> sys          = sss('building');
%> [sysr, V]    = rk(sys,-eigs(sys,4).');
%> [Rv, B_, Sv] = getSylvester(sys, sysr, V);
%
%       One can verify the solution by looking at the norm of the residual:
%> norm(sys.A*V-sys.E*V*Sv-sys.B*Rv)
%
%//Note: RK can return some matrices of the Sylvester equation directly.
% 
% See Also: 
%       rk, porkV, porkW, cure, spark
%
% References:
%       * *[1] Gallivan et al. (2002)*, Sylvester equations and projection
%              based model reduction
%       * *[2] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%       * *[3] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  parsing
if nargin == 3
    type = 'V';
end
if strcmp(type,'W')
    % use the dual system
    sys = sys.'; sysr = sysr.';
end

%%  computations
B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
R = B_\(sys.A*V - sys.E*V*(sysr.E\sysr.A));

if nargout > 2
    S = sysr.E\(sysr.A - sysr.B*R);
end

%% control the accuracy by computing the residual
res = zeros(1,2);
res(1) = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-B_*R);
if nargout > 2
    res(2) = norm(sys.A*V - sys.E*V*S - sys.B*R);
end

%%  Change shape of C_ and Sw 
if strcmp(type,'W'),
    B_ = B_.'; 
    if nargout > 2
        S = S.';
    end
end

%%  Check residuals
if any( res > 1e-6 )
    resMax = max(res);
    if  resMax < 1e-1
        warning('Careful, the problem might be ill conditioned and the results of getSylvester inaccurate (res = %e)',resMax);
    else
        error('The Sylvester equation residual (%e) indicates that getSylvester failed to get the correct solution. Check the condition number of your problem',resMax);
    end
end