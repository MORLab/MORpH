function [p,v] = shamPassive(sys,varargin)
% SHAMPASSIVE - Check if LTI system sys is passive or not by checking if an
%               associated Hamiltonian eigenvalue problem has purely imagina-
%               ry eigenvalues or not. Afterwards, the maximum violation in
%               each band is found by sampling.
%
% Syntax:
%   p = shamPassive(sys)
%   p = shamPassive(sys,opts)
%
% Description:
%       p = shamPassive(sys) constructs a Hamiltonian matrix pencil
%       associated to system sys and checks if purely imaginary eigenvalues
%       are present. Afterwards, the maximum violation in non passive frequency bands are
%       evaluated by sampling.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss object, containing LTI system
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .tol: 	Tolerance for specifiying if the Popov function is
%                       negative or not.
%                       [{1e-8} / positive double]
%           - .tolImag: Tolerance defining if an eigenvalue is tested
%                       for being purely imaginary.
%                       [{1e-1} / positive double]
%
% Output Arguments:
%       - p:    flag specifiying if system is passive or not
%       - v:    Structure containing passivity violation informations
%
% See also:
%       samPassive, HamEig
%
% References:
%       [1] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: 
%           Theory and Applications, 2016.
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
warning('shamPassive uses hamEig, an undocumented mexfile that could change without further notice! Use samPassive for more reliable results.')
[sys,opts] = parseInputs(sys,varargin{:});

% Check if sys is phs object
if isa(sys,'phs') || isa(sys,'phsRed')
    p=true;
    v.violations = [];
    v.frequencies = [];
    v.bands = [];
    v.eigenvectors = [];
    return
end

%% Check D+D'
DDT = sys.D + sys.D.';
eigDDT = eig(DDT);

% Check for positive definiteness
if ~min(eigDDT)>0
    error('shamPassive only works for systems with D+D''>0.')
end

%% Calculate purely imaginary eigenvalues of Hamitlonian matrix pencil
[x_imag,ev_imag,dN0] = hamEig(sys,opts);
v.dN0 = dN0;

%% Check passivity
% Check Popov function at arbitrary frequency
psi = evalG(sys,1i) + evalG(sys,-1i).';
eta = min(real(eig(full(psi))));

p=false;
if isempty(x_imag) && eta > -opts.tol
    p=true;
    v.violations = [];
    v.frequencies = [];
    v.bands = [];
    v.eigenvectors = [];
    return

elseif isempty(x_imag) % No purely imaginary eigenvalues but not passive
    % Define sample points
    samples = logspace(-12,5,500);

    % Initialize minimal eigenvalue vector
    lambda = zeros(length(samples),1);

    % Calculate Popov function in frequency band
    for j = 1:length(samples)
        G=evalG(sys,1i*samples(j));
        Psi=G+G';
        D = eig(full(Psi));
        lambda(j) = min(real(D));
    end

    % Compute minimal eigenvalue
    [m,idx] = min(lambda);
    v.violations = m;
    v.frequencies = samples(idx);
    v.eigenvectors = [];
    v.bands = [0,Inf];
    return
end

%% Check frequency subbands
[f,idx] = sort(x_imag,'ascend');
ev_imag = ev_imag(:,idx);
v.eigenvectors = ev_imag;

% Add zero to frequency vector
if f(1)~=0
    if isrow(f)
        f=f.';
    end
    f=[0;f];
end

% Add new frequency after the last
f(end+1)=10^(2+(floor(log10(f(end)))));

% Midpoint frequencies of each subband
epsilon = zeros(length(f)-1,1);
for i=1:length(epsilon)
    epsilon(i) = (f(i) + f(i+1))/2;
end

% Minimum eigenvalues at midpoint frequencies
eta = zeros(length(epsilon),1);
for i=1:length(eta)
    psi = evalG(sys,1i*epsilon(i)) + evalG(sys,-1i*epsilon(i)).';
    eta(i) = min(real(eig(full(psi))));
end

% Check signs of smallest eigenvalues
signs = sign(eta);
if all(signs>=0)
    p = true;
    v.violations = [];
    v.frequencies = [];
    v.bands = [];
    v.eigenvectors = [];
    return
end

% Calculate subbands
violations = nnz(signs<0);
bands = zeros(violations,2);
j=1;
for i = 1:length(signs)
    if signs(i)<0
        bands(j,1)=f(i);
        bands(j,2)=f(i+1);
        j=j+1;
    end
end

% Merge adjacent passivity violating frequency bands
b=reshape(bands',[],1);
idx = [];
for i=1:length(b)-1
    if b(i)==b(i+1)
        idx = [idx,i,i+1];
    end
end
b(idx)=[];
violations = length(b)/2;
bands = zeros(violations,2);
row = 1;
for i=1:2:numel(b)
    bands(row,1)=b(i);
    bands(row,2)=b(i+1);
    row = row+1;
end

%% Calculate maximum violations
% Initialize
v.violations = zeros(size(bands,1),1);
v.frequencies = zeros(size(bands,1),1);

% Sample over violating frequency bands
for i = 1:size(bands,1)
    % Calculate number of samples
    diff = bands(i,2)-bands(i,1);
    k = floor(log10(diff));
    if k<0
        k=0;
    end
    n = (k+1)*100;

    % Get frequency samples
    if bands(i,1) == 0
        samples = exp(linspace(log(1e-12),log(bands(i,2)),n));
    else
        samples = exp(linspace(log(bands(i,1)),log(bands(i,2)),n));
    end

    % Initialize minimal eigenvalue vector
    lambda = zeros(length(samples),1);

    % Calculate Popov function in frequency band
    for j = 1:length(samples)
        G=evalG(sys,1i*samples(j));
        Psi=G+G';
        D = eig(full(Psi));
        lambda(j) = min(real(D));
    end

    % Compute minimal eigenvalue
    [m,idx] = min(lambda);
    v.violations(i) = m;
    v.frequencies(i) = samples(idx);
end

v.bands = bands;
v.w = x_imag;

%% Calculate slopes
v.slopes = slopes(ev_imag,dN0);

end

%% Auxilliary functions
function [ximag,evimag,dN0] = hamEig(sys,opts)
% Invert D+D'
Dinv = eye(size(sys.D))/(sys.D+sys.D.');

% Construct Hamiltonian system
A=sys.A-sys.B*Dinv*sys.C;
G=-sys.B*Dinv*sys.B.';
Q=sys.C.'*Dinv*sys.C;
N0 = [A,G;Q,-A'];

% Calculate Hamiltonian eigenvalues
x=HamEig(full(A),full(G),full(Q));

% Remove almost zero real and imaginary parts
for i=1:length(x)
    if abs(real(x(i)))<1e-6
        x(i) = 0+1i*imag(x(i));
    end
end

for i=1:length(x)
    if abs(imag(x(i)))<1e-6
        x(i) = real(x(i));
    end
end

% Check for duplicate eigenvalues
idx_ev=[];
for i=1:length(x)
    active = x(i);
    for j=i+1:length(x)
        if abs(active-x(j))<1e-6
            x(j)=NaN;
            idx_ev(end+1)=j;
        end
    end
end

x=x(~isnan(x));

% Get derivative of Ah
Dinv2 = eye(size(sys.D))/Dinv;
dA = -2*sys.B*Dinv2*sys.C;
dG = -2*sys.B*Dinv2*sys.B.';
dQ = 2*sys.C.'*Dinv2*sys.C;
dN0 = [dA,dG;dQ,-dA.'];

% Get eigenvalues near the imaginary axis
ximag = x(abs(real(x))<opts.tolImag);

% Sort them such that mirror images are stacked together
ximag = sort(ximag);

% Filter purely imaginary eigenvalues
ximag = pureImagEig(ximag);

evimag = zeros(size(N0,1),length(ximag));
for i = 1:length(ximag)
    vec = null(full(N0)-ximag(i)*eye(size(N0)));
    evimag(:,i) = vec(:,1)/norm(vec(:,1));
end

% Return only imaginary part
ximag = imag(ximag);
end


function [ximag] = pureImagEig(x)
% Create empty arrays and take only eigenvalues into account with
% positive imaginary part
ximag = [];
x = x(imag(x)>=0);

% Iterate over all eigenvalues
i=1;
while i<=length(x)
    % Get active eigenvalue
    active_val = x(i);

    % Get imaginary and real part
    real_x = real(active_val);
    imag_x = imag(active_val);

    % Check if mirror image is present
    mirror_image_present = false;
    for j=i+1:length(x)
        if abs(-real_x-real(x(j))) < 1e-10
            if abs(imag_x-imag(x(j))) < 1e-10
                mirror_image_present = true;
                i=i+1;
                break;
            end
        end
    end

    % If not present, eigenvalue is purely imaginary
    if ~mirror_image_present
        ximag(end+1) = active_val;
    end

    % Raise counter
    i=i+1;
end
end

function s = slopes(Xr,dN0)
n = size(dN0,1)/2;
J = [zeros(n),eye(n);-eye(n),zeros(n)];

s = zeros(size(Xr,2),1);
for i = 1:length(s)
    s(i) = real((1i*Xr(:,i)'*J*Xr(:,i))/(Xr(:,i)'*J*dN0*Xr(:,i)));
end
end

function G_si = evalG(sys,s_i)
% Evaluates the transfer function G(s) at certain shifts s_i.
res = (s_i*sys.E-sys.A)\sys.B;
G_si = sys.C*res+sys.D;
end

function [sys,opts] = parseInputs(sys,varargin)
% Check if opts is provided
if nargin < 2
    opts = struct();
else
    opts = varargin{1};
end

% Set default values
OptsAdmissible.tol = 1e-8;
OptsAdmissible.tolImag = 1e-1;

% Set options for sadpa/samdp
OptsAdmissible.eig.nwanted = 50;       % Number of wanted poles
OptsAdmissible.eig.strategy = 'LR';    % dominance measure according to [1]
OptsAdmissible.eig.displ = 0;          % dont show convergence results
OptsAdmissible.eig.use_lu = 0;

opts = phsMOR_parseOpts(opts, OptsAdmissible);

if ~isstable(ss(sys))
    error('System is not asymptotically stable.')
end

% Check if sys is dae
if svds(sparse(sys.A),1,'smallest') < 1e-12
    error('shamPassive does not support DAE systems.')
end

if ~isempty(sys.E)
    % Make system explicit
    sys.A = sys.E\sys.A;
    sys.B = sys.E\sys.B;
    sys.E = eye(size(sys.A));
else
    sys.E = eye(size(sys.A));
end

end