function sysr = dszm(sys, r, varargin)
% dszm - Obtain a passive reduced order model by rational
%               interpolation at spectral zeros.
%
% Syntax:
%   sysr = dszm(sys, r)
%   sysr = dszm(sys, r, Opts)
%
% Description:
%       sysr = dszm(sys, r, Opts) reduces the passive system sys by
%       utilizing the interpolatory Arnoldi method [4] with a set of dominant
%       spectral zeros and the according residues. The arising Hamiltonian
%       eigenvalue problem is solved by SADPA/SAMDP [5,6].
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss object, containing passive LTI system
%       - r:        desired reduced order
%
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .tol:                 If D=0 it will be subsituted by tol*eye(m).
%                                   [{1e-12} / positive double]
%           - .explicit:            Define if explicit system should be
%                                   returned.
%                                   [{true} / false]
%           - .makePH:              Make reduced system port-Hamiltonian.
%                                   [{false} / true]
%           - .checkPassivity:      Check passivity of sys.
%                                   [{false} / true]
%           - .projectiveMOR.*: 	Other options that will be passed on to the
%                                   used method (sssMOR->projectiveMOR);
%                                   Please refer to documentation of the respective
%                                   algorithm.
%           - .arnoldi.*:           Other options that will be passed on to the
%                                   used method (sssMOR->arnoldi);
%                                   Please refer to documentation of the respective
%                                   algorithm.
%           - .samPassive.*:        Other options that will be passed on to the
%                                   used method (Utility/passivity_check/samPassive);
%                                   Please refer to documentation of the respective
%                                   algorithm.
%           - .ss2phs.*:            Other options that will be passed on to the
%                                   used method (Utility/ss2phs);
%                                   Please refer to documentation of the respective
%                                   algorithm.
%           - .sadpa.*:             Other options that will be passed on to the
%                                   third-party software 'sadpa';
%                                   Pleae refer to documentation of the respective
%                                   algorithm.
%           - .samdp.*:             Other options that will be passed on to the
%                                   third-party software 'samdp';
%                                   Pleae refer to documentation of the respective
%                                   algorithm.
%
% Output Arguments:
%       - sysr: 	ssRed (phsRed) object, containing passive reduced LTI (port-Hamiltonian) system
%
% See Also:
%       arnoldi, projectiveMOR, sadpa, samdp, demo_dszm
%
% References:
%       [1] R. Ionutiu, J. Rommes, and A. C. Antoulas. " Passivity-Preserving 
%           Model Reduction Using Dominant Spectral-Zero Interpolation." 
%           In: IEEE Transactions on Computer-Aided Design of Integrated 
%           Circuits and Systems 27.12 (2008), pp. 2250–2263.
%       [2] D.C. Sorensen. Passivity preserving model reduction via interpolation 
%           of spectral zeros. In 2003 European Control Conference (ECC), 
%           pp. 974–978, 2003.
%       [3] A. C. Antoulas. A new result on passivity preserving model reduction.
%           Systems & Control Letters, 54:361–374, 2005.
%       [5] J. Rommes and N. Martins. “Efficient Computation of
%           Multivariable Transfer Function Dominant Poles Using
%           Subspace Acceleration.” In: IEEE Transactions on Power
%           Systems 21.4 (2006), pp. 1471–1483.
%       [6] J. Rommes and N. Martins. “Efficient computation of
%           transfer function dominant poles using subspace acceleration.” 
%           In: IEEE Transactions on Power Systems 21.3 (2006), pp. 1218–1226.
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Maximilian Bonauer, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%               This function uses the third-party functions 'sadpa and 'samdp' [5,6], 
%               please refer to their specific license.
%-----------------------------------------------------------------------

%% Input parsing
[r,m,D,posDefD,Opts] = parseInputs(sys, r, varargin{:});

%% Compute associated Hamiltonian system
% Inverse D+D'
Delta = sparse((D+D')\eye(size(D)));

% Set up Hamiltonian system
if posDefD
    Ah = [sys.A-sys.B*Delta*sys.C, -sys.B*Delta*sys.B'; sys.C'*Delta*sys.C, -sys.A'+sys.C'*Delta*sys.B'];
    Eh = blkdiag(sys.E,sys.E');
    Bh = [sys.B;sys.C']*Delta;
    Ch = -Delta*[sys.C, sys.B'];
else
    Ah = [blkdiag(sys.A, -sys.A'), [sys.B; -sys.C'];...
        [sys.C, sys.B'], D+D'];
    Eh = blkdiag(sys.E, sys.E', zeros(m));
    Bh = sparse([sys.B; -sys.C'; zeros(m)]*Delta);
    Ch = sparse(-Delta*[sys.C, sys.B', zeros(m)]);
end

%% Compute spectral zeros
if m==1
    spectralZeros = sadpa(Ah, Eh, Bh, Ch', Delta, 1i, Opts.sadpa);
else
    spectralZeros = samdp(Ah, Eh, Bh, Ch', Delta, 1i, Opts.samdp);
end

%% Filter spectral zeros
% Get all spectral zeros with negative real part
negativeSpectralZeros = spectralZeros(spectralZeros<0);

% Remove duplicate entries
negativeSpectralZeros = unique(negativeSpectralZeros,'stable');

% Get the most dominant negative spectral zeros
shifts = [];
i=1;
while i<length(negativeSpectralZeros)
    if isreal(negativeSpectralZeros(i))
        shifts(i) = negativeSpectralZeros(i);
        i=i+1;
    elseif abs(imag(negativeSpectralZeros(i))) < 1e-10
        shifts(i) = real(negativeSpectralZeros(i));
        i=i+1;
    else
        shifts(i) = negativeSpectralZeros(i);
        shifts(i+1) = conj(negativeSpectralZeros(i));
        i=i+2;
    end

    if length(shifts)>=r
        break
    end
end

% Check if any spectral zeros were found
if length(shifts)<=0
    error('No spectral zeros for interpolation were found.')
end

% Check if enough spectral zeros were found
if length(negativeSpectralZeros) < r
    r_old = r;
    r = length(negativeSpectralZeros);
    if m==1
        warning('Reduced order was decreased from %i to %i because the number of negative spectral zeros found by sadpa was too small.',r_old,r)
    else
        warning('Reduced order was decreased from %i to %i because the number of negative spectral zeros found by samdp was too small.',r_old,r)
    end
end

% Check if the reduced order increased
if length(shifts) > r
    r = length(shifts);
    warning('Reduced order increased by %i to %i due to distribution of spectral zeros.',m,r)
end

%% Projection
% Get projection matrices V and W
V = arnoldi(sys.E, sys.A, sys.B, shifts, Opts.arnoldi);
W = arnoldi(sys.E.', sys.A.', sys.C.', -conj(shifts), Opts.arnoldi);

% Ensure E=I
if Opts.explicit
    [L,U] = lu(W'*sys.E*V);
    V = V/U;
    W = (L\W')';
end

% Compute reduced order model
sysr = projectiveMor(sys, V, W, Opts.projectiveMOR);

% Stability check
if ~isstable(ss(sysr))
    error('Reduced order model is unstable.')
end

% Passivity check
passive = samPassive(sysr,Opts.samPassive);
if ~passive
    error('Reduced order model is not passive.')
end

% Transform to pH form
if Opts.makePH
    sysr = ss2phs(sysr, Opts.ss2phs);
    sysr.method = @dszm;
    sysr.parameters = Opts;
    sysr.parameters.spectralZeros = [shifts;-conj(shifts)];
else
    % Create ssRed object
    Opts.spectralZeros = [shifts;-conj(shifts)];
    sysr = ssRed(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,'dszm',Opts,sys);
end

end

%% Auxilliary function
function [r,m,D,posDefD,Opts] = parseInputs(sys, r, varargin)

% Check if Opts is provided
if nargin < 2
    error('Not enough input arguments.')
elseif nargin < 3
    Opts = struct();
else
    Opts = varargin{1};
end

% Check input arguments
isPHS = false;

% Check if FOM is pH
if isa(sys,'phs') || isa(sys,'phsRed')
    sys = phs2sss(sys);
    isPHS = true;
end

% Check if sys is dae
if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
    error('dszm currently only supports ODE systems.')
end

% Number of inputs/outputs
m = size(sys.B,2);

% Check if SISO or MIMO
if m>1
    warning('Dominant spectral zero method only works reliable for SISO systems.')
end

% Check if r is multiple of m
if mod(r,m) ~= 0
    error('Reduced order r must be a multiple of number of inputs m.')
end

% Set default values and asmissible value sets
OptsAdmissible.tol = 1e-12;
OptsAdmissible.explicit = {true,false};
OptsAdmissible.checkPassivity = {false,true};
OptsAdmissible.makePH = {false,true};
OptsAdmissible.projectiveMOR = struct();
OptsAdmissible.arnoldi = struct();
OptsAdmissible.samPassive.plot = false;
OptsAdmissible.ss2phs = struct();

% Set default options for sadpa/samdp
dszOpts.nwanted = 5*r;      % Compute 5*r spectral zeros because sometimes sadpa/samdp return not enough spectral zeros
dszOpts.strategy = 'LR';    % dominance measure according to [1]
dszOpts.displ = 0;          % dont show convergence results

if m == 1
    OptsAdmissible.sadpa = dszOpts;
else
    OptsAdmissible.samdp = dszOpts;
end

Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

% Check for passive system
if ~isPHS && Opts.checkPassivity
    passive = samPassive(sys,Opts.samPassive);
else
    passive = true;
end

if ~passive
    error('Dominant spectral zero method is only applicable to passive systems.')
end

% Opts.explicit must be true if .makePH is true
if Opts.makePH && ~Opts.explicit
    Opts.explicit = true;
    warning('ss2phs needs an explicit reduced order model. Setting Opts.explicit to true.')
end

% Check for D+D' > 0 and correct if necessary
ev = eig(sys.D+sys.D');
if all(ev>0)
    D = sys.D;
    posDefD = true;
else
    D = sys.D + Opts.tol*eye(size(sys.D));
    posDefD = false;
end

% Check third-party software
if m==1
    thirdPartyCheck('SADPA');
else
    thirdPartyCheck('SAMDP');
end

end