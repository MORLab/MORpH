function hasStaircase = staircaseCheck(sys, varargin)
% STAIRCASECHECK - Checks if a pHDAE is in staircase form
%
% Syntax:
%   hasStaircase = staircaseCheck(sys)
%   hasStaircase = staircaseCheck(sys, Opts)
%
% Description:
%                
% The pH system
%     E*dx/dt = (J-R)*x(t)  + (G-P)*u(t),
%           y = (G+P)'*x(t) + (S+N)*u(t),
%
% with partitioned state vector x = [x1;x2;x3;x4] in R^(n1+n2+n3+n4) is in staircase form, if
%
% (i)   The system satisfies the pH structural constraints.
% (ii)  The system matrices have the following zero-patterns:
%      [E11  0  0 0]      [J11 J12 J13 J14]      [R11 R12 R13  0 ]      [G1]      [P1]
%  E = [ 0  E22 0 0], J = [J21 J22 J23  0 ], R = [R21 R22 R23  0 ], G = [G2], P = [P2]
%      [ 0   0  0 0]      [J31 J32 J33  0 ]      [R31 R32 R33  0 ]      [G3]      [P3]
%      [ 0   0  0 0]      [J41  0   0   0 ]      [ 0   0   0   0 ]      [G4]      [0 ]
%
% (iii) rank(E11) = n1, rank(E22) = n2
% (iv)  rank(J14) = n1 = n4
% (v)   rank(J33-R33) = n3
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object 
%       *Optional Input Arguments:*
%       - Opts:
%             - rankTol: Tolerance for rank decisions
%                        [{1e-12} / positive double]
%             - verbose: If true, the function throws an error with more details
%                        if the model is not in staircase form; otherwise the 
%                        output is simply set to false
%                        [{false} / true]
%
% Output Arguments:
%       - hasStaircase: 'true' if system has staircase form, 'false' otherwise
%
% See Also:
%       staircaseDims, toStaircase
%
% References:
%       [1] F. Achleitner, A. Arnold, and V. Mehrmann. Hypocoercivity and controllability in linear
%           semi-dissipative ODEs and DAEs. ZAMM Z. Angew. Math. Mech., 2021.
%
%-----------------------------------------------------------------------
% This file is part of 
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a> 
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Input parsing
narginchk(1,2);
[sys, Opts] = parseInputs(sys, varargin{:});

%% Check staircase conditions
hasStaircase = true;
try

    errMsg = 'The provided system does not have staircase form.\n';

    % Check if system is in scaled energy coordinates
    assert(isequal(sys.Q,eye(sys.dim)),'phs:staircaseCheck:notInScaledEnergyCoordinates',...
        strcat(errMsg,'System is not in scaled energy coordinates (Q~=I).\n'));

    % Get dimensions (assuming theat sys has staircase form)
    dims = staircaseDims(sys);
    n1 = dims(1); n2 = dims(2); n3 = dims(3); n4 = dims(4);

    assert(all(dims >= 0) && n1 == n4,'phs:staircaseCheck:stateDecompositionFailed',...
        strcat(errMsg,['Decomposition of the state vector failed.\n' ...
        'This is an indication that the system does not have staircase form.\n']));

    % Check zero patterns
    errMsgZeroPat = 'Matrix %s does not have the required zero pattern.\n';

    if n1 + n2 > 0
        assert(isequal(sys.E-blkdiag(sys.E(1:n1,1:n1),sys.E(n1+1:n1+n2,n1+1:n1+n2),zeros(n3+n4)), zeros(size(sys.E))),...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgZeroPat,'E')));
    end

    if n1 > 0
        assert(isequal(sys.J(end-n4+1:end,n1+1:end), zeros(n4,n2+n3+n4)) ...
            && isequal(sys.J(n1+1:end,end-n4+1:end), zeros(n2+n3+n4,n4)), ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgZeroPat,'J')));
        assert(isequal(sys.R-blkdiag(sys.R(1:n1+n2+n3,1:n1+n2+n3),zeros(n1)), zeros(size(sys.R))), ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgZeroPat,'R')));
        assert(isequal(sys.P(end-n4+1:end,:), zeros(n4,size(sys.P,2))),...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgZeroPat,'P')));
    end

    % Check rank conditions
    errMsgRankCond = ['Matrix %s does not have full rank. \n' ...
        'Consider changing the tolerance for the smallest singular value (currently %.2d).\n'];

    if n1 > 0
        assert(svds(sys.E(1:n1,1:n1), 1, 'smallest') > Opts.rankTol, ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgRankCond,'E11',Opts.rankTol)));
        assert(svds(sys.J(end-n4+1:end,1:end), 1, 'smallest') > Opts.rankTol, ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgRankCond,'J41',Opts.rankTol)));
    end

    if n2 > 0
        assert(svds(sys.E(n1+1:n1+n2,n1+1:n1+n2), 1, 'smallest') > Opts.rankTol, ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgRankCond,'E22',Opts.rankTol)));
    end

    if n3 > 0
        assert(svds(sys.J(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3)-sys.R(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3), 1, 'smallest') > Opts.rankTol, ...
            'phs:staircaseCheck:wrongZeroPattern',...
            strcat(errMsg,sprintf(errMsgRankCond,'J33-R33',Opts.rankTol)));
    end
    
catch ME
    if Opts.verbose
        rethrow(ME)
    else
        hasStaircase = false; 
    end
end

end

function [sys, Opts] = parseInputs(sys, varargin)
    % Check phs input type
    if ~isa(sys,'phs')
        error('phs:staircaseCheck:wrongInput', 'Model is not an object of the phs-class.');
    end
    
    % Opts
    if ~isempty(varargin) && isstruct(varargin{end})
        Opts = varargin{end};
    else
        Opts = struct();
    end
    
    % Option parsing
    OptsAdmissible.rankTol = 1e-12;
    OptsAdmissible.verbose = false;
    
    Opts = phsMOR_parseOpts(Opts,OptsAdmissible);
end