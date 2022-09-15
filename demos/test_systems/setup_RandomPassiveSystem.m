function [sys,sysPH] = setup_RandomPassiveSystem(n,m,Opts)
% SETUP_RANDOMPASSIVESYSTEM - Construct a passive system of order n with in- 
% and output dimension m
%
% Syntax:
%   [sys,sysPH] = setup_RandomPassiveSystem(n,m,Opts)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - n:	system dimension
%       - m:    in- and output dimension of the system
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .rankW:	Rank of W or R. Default: n-m
%           - .ft:      Define if feed-through is present or not
%                       Default: true
%           - .P:       Define if P is random or zero.
%                       Default: 'zero'
%           - .phs.*:   Other options that will be passed on to the
%                       used method (@phs/phs);
%                       Refer to documentation of the respective
%                       algorithm.
% 
% Output Arguments:
%       - sys:      sss-object
%       - sysPH: 	phs-object
%
% Examples:
%       The following code generates a passive SISO system of order 10
%
%       sys = setup_RandomPassiveSystem(10,1)
%
% See Also:
%       phs, setup_RandomStaircasePHDAE
%
%-----------------------------------------------------------------------
% This file is part of 
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Maximilian Bonauer
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a> 
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Input parsing
if nargin < 2 
    error('Not enough input arguments.')
end

if nargin < 3
    Opts = struct();
end

if ~isnumeric(n)
    error('Input dimension n is not numeric.')
end

if ~isnumeric(m)
    error('Input dimension m is not numeric.')
end

OptsAdmissible.phs.inputValidation = {true,false};
OptsAdmissible.phs.inputTolerance = 1e-8;
OptsAdmissible.rankW = -1;
OptsAdmissible.ft = {true,false};
OptsAdmissible.P = {'rand','zero'};

Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% 
if Opts.rankW == -1
    Opts.rankW = n;
end

if Opts.rankW < 1
    error('Rank of W must be greater than 0.')
end

if Opts.rankW > n
    error('Rank of W can not be greater than %2d.',n)
end

if strcmp(Opts.P,'rand') && ~Opts.ft
    warning('P can only be nonzero if a feed-through is present.')
    Opts.P = 'zero';
end

%% Dimensions
% Dimension of matrix Gamma
dims.nJ = n*(n-1)/2;
dims.nG = n*m;
dims.nN = m*(m-1)/2;
dims.nG = dims.nJ + dims.nG + dims.nN;

%% Get parameter vectors
thetaG = rand(dims.nG,1);
if strcmp(Opts.P,'rand')
    thetaW = rand(n+m,Opts.rankW);
else
    thetaR = rand(n,Opts.rankW);
    thetaS = rand(m,m);
end

%% Create random pH system
alpha=0.05; % Tuneable parameter
if strcmp(Opts.P,'rand')
    W = alpha*(thetaW*thetaW');
end
Gamma = vtsu(thetaG)'-vtsu(thetaG);

J = Gamma(1:n,1:n);

if strcmp(Opts.P,'rand')
    R = W(1:n,1:n);
else
    R = alpha*(thetaR*thetaR');
end

Q = randomPD(n,rand(1));
G = Gamma(1:n,n+1:end);
E = eye(n);

if Opts.ft
    if strcmp(Opts.P,'rand')
        P = W(1:n,n+1:end);
        S = W(n+1:end,n+1:end);
    else
        P = zeros(n,m);
        S = thetaS*thetaS';
    end
    
    N = Gamma(n+1:end,n+1:end);
else
    P = zeros(n,m);
    S = zeros(m,m);
    N = zeros(m,m);
end

if n <= 0
    sysPH = phs(0,0,0,zeros(1,m),1,zeros(1,m),S,N);
else
    sysPH = phs(J,R,Q,G,E,P,S,N,Opts.phs);
end

sys = sss(sysPH);

end

function A = randomPD(n,eigen_mean)
Q = randn(n,n);
Q=orth(Q);
% Q = rand(n,n);
% can be made anything, even zero 
% used to shift the mode of the distribution
A = Q'*diag(abs(eigen_mean+randn(n,1)))*Q;
end