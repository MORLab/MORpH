function [sys,sysPH,cur_vio] = setup_RandomSystem(n,m,des_vio,Opts)
% SETUP_RANDOMSYSTEM - Construct a non-passive system of order n with in- 
% and output dimension m and passivity violation des_vio
%
% Syntax:
%   sys = setup_RandomSystem(n,m,des_vio)
%   [sys,sysPH,cur_vio] = setup_RandomSystem(n,m,des_vio,Opts)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - n:        system dimension
%       - m:        in- and output dimension of the system
%       - des_vio:  desired value of the passivity violation
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .filter: 	        Filter width defining the width of the
%                               passivity violation. Default: 0
%           - .rankW: 	        Rank of W or R. Default: n-m
%           - .ft:              Define if feed-through is present or not
%                               [{true} / false]
%           - .P:               Define if P is random or zero.
%                               [{'zero'} / 'rand']
%           - .samPassive.*:    Other options that will be passed on to the
%                               method (Utility/passivity_check/samPassive);
%                               Please refer to documentation of the respective
%                               algorithm.
%           - .phs.*:           Other options that will be passed on to the
%                               method (@phs/phs);
%                               Please refer to documentation of the respective
%                               algorithm.
% 
% Output Arguments:
%       - sys:      sss-object of non-passive system
%       - sysPH: 	phs-object of random pH system used for initialization
%       - cur_vio: 	actual value of passivity violation
%
% Examples:
%       The following code computes a non-passive SISO system of order 10
%       with a maximum passivity violation of -1e-1.
%
%       sys = setup_RandomSystem(10,1,-1e-1)
%
% See Also:
%       samPassive, setup_RandomPassiveSystem
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
    des_vio = -1e-6;
end

if nargin < 4
    Opts = struct();
end

if ~isnumeric(n)
    error('Input dimension n is not numeric.')
end

if ~isnumeric(m)
    error('Input dimension m is not numeric.')
end

if ~isnumeric(des_vio)
    error('Value for violation is not numeric.')
end

if des_vio >= 0
    error('Passivity violation must be less than zero.')
end

OptsAdmissible.phs.inputValidation = {true,false};
OptsAdmissible.phs.inputTolerance = 1e-8;
OptsAdmissible.samPassive.plot = false;
OptsAdmissible.samPassive.outputVolume = 'worst';
OptsAdmissible.filter = {0,1,2,3,4,5,6,7,8,9,10,11,12};
OptsAdmissible.rankW = -1;
OptsAdmissible.ft = {true,false};
OptsAdmissible.P = {'zero','rand'};

Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% Not optional values
tol = 10^(floor(real(log10(des_vio)))-3);
maxit = 2000;

% 
if Opts.rankW == -1
    Opts.rankW = n-m;
end

if Opts.rankW < 1
    error('Rank of W must be greater than 0.')
end

if Opts.rankW > n-m
    error('Rank of W can not be greater than %2d.',n-m)
end

if strcmp(Opts.P,'rand') && ~Opts.ft
    warning('P can only be nonzero if a feed-through is present.')
    Opts.P = 'rand';
end

%% Dimensions
% Dimension n without filter attachment
n=n-2*m;

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
    thetaS = rand(m,m)-des_vio;
end

%% Create random pH system
if strcmp(Opts.P,'rand')
    W = thetaW*thetaW';
end
Gamma = vtsu(thetaG)'-vtsu(thetaG);

J = Gamma(1:n,1:n);

if strcmp(Opts.P,'rand')
    R = W(1:n,1:n);
else
    R = thetaR*thetaR';
end

Q = randomPD(n,rand(1));
G = Gamma(1:n,n+1:end);
Z = randomPD(n,rand(1));
E = (Z/Q)';

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
sysp=phs2sss(sysPH);

%% Construct filter and define sampling band
MSD = ss(setup_MassSpringDamperSystem(2,1,10^-Opts.filter,1,'SISO'));
M = MSD;
for i=1:m-1
    M = blkdiag(M,MSD);
end
filter = 0.5*sss(ss(M));
w = logspace(-Opts.filter-2,2,1000);
[~,v]=samPassive(sysp,w,Opts.samPassive);
d0 = v.violations;

%% Calculate initial guess for x        
x = d0-des_vio;
xlb = 0;
xub = max([20,3*x]);
x = 0.5*(xlb+xub);
sys = sysp-x*filter;

%% Make sys non-passive
% Calculate current violation
[~,v]=samPassive(sys,w,Opts.samPassive);
cur_vio = v.violations;

i = 1;
while abs(cur_vio-des_vio) > tol && i < maxit
    if i ~= 1
        % Compute new x and corresponding sys
        x = 0.5*(xlb+xub);
        sys = sysp-x*filter;
    end
    
    % Get passivity violation
    [~,v]=samPassive(sys,w,Opts.samPassive);
    cur_vio = v.violations;
    
    % Update boundaries
    if cur_vio > des_vio
        xlb = x;
    else
        xub = x;
    end
    i=i+1;
end

% Check if algorithm was successful
if abs(cur_vio-des_vio) > tol && i == maxit
    warning('Algorithm did not converge. Resulting passivity violation is %2.2f instead of %2.6f',cur_vio,des_vio)
end

end

function A = randomPD(n,eigen_mean)
Q = randn(n,n);
Q=orth(Q);
% Q = rand(n,n);
% can be made anything, even zero 
% used to shift the mode of the distribution
A = Q'*diag(abs(eigen_mean+randn(n,1)))*Q;
end