function [redSys, residual, shifts, Dr] = trksmPH(sys,  rMax, varargin)
% TRKSMPH:  Adaptive Tangential Rational Krylov Subspace Method for Port-Hamiltonian Systems
%
% Syntax:
%       redSys         = trksmPH(sys, rMax)
%       redSys         = trksmPH(sys, rMax, Opts)
%       redSys         = trksmPH(sys, rMax, sMin, sMax)
%       redSys         = trksmPH(sys, rMax, sMin, sMax, Opts)
%       redSys         = trksmPH(sys, rMax, sMin, sMax, lFix)
%       redSys         = trksmPH(sys, rMax, sMin, sMax, lFix, Opts)
%
% Description:
%
% TRKSMPH iteratively builds up a reduced-order pHODE model. In each
% iteration, a new interpolation point is added in areas, where the
% approximation quality (measured by a residual term) is still poor. This
% way, the reduced dimension increases in each iteration until the process
% converges.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:	    phs object of the large-scale model
%       - rMax:     maximum size of the reduced system
%
%       *Optional Input Arguments:*
%       - sMin, sMax:   Initial shifts
%                       [{sMin=1e-5,sMax=1e5} / complex double with real(sMin) < real(sMax)]
%       - lFix:         Use 'lFix' tangential directions for each shift; for lFix == 0,
%                       the number of used tangential directions will be chosen automatically
%                       [{0} / positive integer < number of inputs]
%       - Opts:     Structure with execution parameters
%           - .tol:              relative H2-norm between two iteration steps. 
%                                If norm(sysnew-sysold)/norm(sysnew) < tol,
%                                the process will be ended
%                                [{1e-3} / positive double]
%           - .minSamplePoints:  Minimum number of sampling points used for
%                                detecting new shifts (default: 500)
%                                [{500} / positive integer]
%           - .isReal:           variable, that defines whether only real or 
%                                also complex poles are allowed. 
%                                1: only real poles 
%                                0: both real and complex poles
%                                [{0} / 1]
%           - .border:           variable, that defines whether poles are chosen 
%                                on convex hull or over convex set
%                                1: only on hull, 0: also within set. 
%                                [{1} / 0]
%           - .log:              only used if border = 0. Defines sampling method over convex set.
%                                1: logarithmic sampling, 0: linear sampling
%                                [{1} / 0]
%
% Output Arguments:
%       - redSys:       reduced pH system (phsRed object)
%       - residual:     residual for each step
%       - shifts:       interpolation points
%       - Dr:           tangential directions
%
% See Also:
%       irkaPH, arnoldiPH
%
% References:
%       [1] V. Druskin and V. Simoncini. “Adaptive rational Krylov
%           subspaces for large-scale dynamical systems.” In: Systems
%           & Control Letters 60.8 (2011), pp. 546–560.
%       [2] V. Druskin, V. Simoncini, and M. Zaslavsky. “Adaptive
%           Tangential Interpolation in Rational Krylov Subspaces for
%           MIMO Dynamical Systems.” In: SIAM Journal on Matrix
%           Analysis and Applications 35.2 (2014), pp. 476–498.
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Nora Reinbold, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Input Parsing

narginchk(2,6);
[rMax, sMin, sMax, lFix, Opts] = parseInputs(sys, rMax, varargin{:});

% Bring to scaled energy coordinates
sys = scaling(sys);

% Original matrices
A = sys.J - sys.R;
B = sys.G - sys.P;
E = sys.E;

%% Compute initial Eta

% Include complex conjugate shifts
if Opts.isReal
    Eta = real([sMin, sMax]);
else
    Eta = unique(sort([[sMin, sMax],[sMin, sMax]'.'], 'ComparisonMethod', 'real'));
end

%% Compute initial ROM with block interpolation at sMin
shifts = sMin;
Dr = eye(size(B,2));
Dr_out{1} = Dr;
Vr_tilde = (A-sMin*E)\B;
if imag(sMin)
    Dr = [Dr , Dr];
    shifts(2) = sMin';
    Dr_out{2} = conj(Dr);
    Vreal = real(Vr_tilde);
    Vimag = imag(Vr_tilde);
    Vr_tilde = [Vreal, Vimag];
end

% Orthogonalize V2
[V2Orth, Rm] = qr(Vr_tilde, 0); %QR-Decomposition for orthogonalization purposes

% Compute ROM
redOrder = size(V2Orth, 2);
[Jr, Rr, Qr, Gr, Er, Pr, Sr, Nr] = structurePreservation(sys, V2Orth, Opts.structurePreservation);
Ar = Jr - Rr;
Br = Gr - Pr;

% Create strictly-proper subsystem
syssprOld = dss(full(Ar), full(Br), full((Gr + Pr)'), 0, full(Er));

%initialize residual
residual = [];

while redOrder < rMax

    % Compute eigenvalues of proper part of ROM
    Lambda = eig(full(Ar), full(Er));
    EtaSet = sort([-Lambda.', Eta], 'ComparisonMethod', 'real');

    % Compute convex hull of EtaSet
    if any(imag(EtaSet)) && (rMax-redOrder) > 1
        index = convhull(real(EtaSet), imag(EtaSet));
        hull = sort(EtaSet(index), 'ComparisonMethod', 'real'); %border of convex hull
    else
        hull = [EtaSet(1),EtaSet(end)];
    end

    % Analyze only points with imag() >= 0 (residual is symmetric w.r.t. real axis)
    hull = hull(imag(hull) >= 0);

    % Create sampling points on the convex hull of EtaSet
    sBorder = sampleBorder(hull, Opts);

    % Include sampling points of the interior if desired
    sInner = [];
    if Opts.border == 0 && any(imag(hull))
        sInner = sampleInterior(hull, Opts);
    end

    samples = sort([sBorder, sInner], 'ComparisonMethod', 'real');

    % Compute new shift based on residual evaluations at sampling points
    [sNew, resNew] = newShift(Ar, B, Br, E, Er, Dr, Rm, samples, V2Orth);

    shifts = [shifts, sNew];

    residual = [residual, resNew];

    % Compute tangential directions in MIMO case
    if sys.isMIMO == 1

        % Computation of the right singular vectors d
        [~, L]=qr(B-E*V2Orth*(Er\(Br)), 0);
        [~, singval, dr]= svd(L*((Dr/Rm)*((Ar-sNew*Er)\(Br))-eye(size(B,2))));

        % Discard non-dominant directions
        if lFix > 0 % choose fixed number of tangential directions
            l = lFix;
        else % Choose number of tangential directions automatically
            sval = diag(singval);
            svalDom = sval(sval > (1/10)*sval(1));
            l = length(svalDom);
        end

        % Make sure that there is enough space
        if imag(sNew)
            lMax = floor((rMax-redOrder)/2);
        else
            lMax = floor((rMax-redOrder));
        end
        dr = dr(:,1:min([l,lMax]));
    else
        dr = 1;
    end

    Dr_out{end+1} = dr;

    % Compute V
    if sum(shifts == sNew) >= 2     % If shift has already been used, compute next Krylov-Direction
        mult = sum(shifts == sNew);
        VrNew = (((A-sNew*E)^mult)\B)*dr;
    else
        VrNew = ((A-sNew*E)\B)*dr;
    end

    % if shift is complex, include complex conjugate shift and split V and dr
    if imag(sNew)
        Vreal = real(VrNew);
        Vimag = imag(VrNew);
        VrNew = [Vreal, Vimag];
        Dr = [Dr, real(dr), imag(dr)];
        shifts(end+1) = sNew';
        Dr_out{end+1} = conj(dr);
    else
        dr = real(dr);
        Dr = [Dr, dr];
    end
    Vr_tilde = [Vr_tilde, VrNew];

    % Orthogonalize V2
    [V2Orth, Rm] = qr(Vr_tilde, 0); % QR-Decomposition for orthogonalization purposes

    % Compute ROM
    redOrder = size(V2Orth,2);
    [Jr, Rr, Qr, Gr, Er, Pr, Sr, Nr] = structurePreservation(sys, V2Orth, Opts.structurePreservation);
    Ar = Jr - Rr;
    Br = Gr - Pr;

    % Create strictly-proper subsystem
    syssprNew = dss(full(Ar), full(Br), full((Gr + Pr)'), 0,full(Er));

    % Compute relative error and compare with tolerance
    err = norm(syssprNew - syssprOld)/norm(syssprNew);
    if err < Opts.tol
        break;
    end

    syssprOld = syssprNew;

    %     % Plot for debugging purposes
    %     SmBorder = [sBorder, (sBorder').'];
    %     SmInner = [sInner, (sInner').'];
    %     figure;
    %     semilogx(real(SmBorder), imag(SmBorder),'r*'); hold on;
    %     semilogx(real(SmInner), imag(SmInner),'b*')
    %     title('Mirrored spectral region')
    %     xlabel('Re(S_m)')
    %     ylabel('Im(S_m)')

end

%% Create phsRed object
redSys = phsRed(Jr, Rr, Qr, Gr, Er, Pr, Sr, Nr);

redSys.method = @trksmPH;
redSys.parameters = Opts;
redSys.parameters.sMin = sMin;
redSys.parameters.sMax = sMax;
redSys.parameters.lFix = lFix;
redSys.parameters.rMax = rMax;

redSys.info = struct();
redSys.info.residual = residual;
redSys.info.shifts = shifts;
redSys.info.Dr = Dr_out;

end

%% ======================== AUXILIARY FUNCTIONS ===========================

function [shift, residual] = newShift(Apr, Bp, Bpr, Ep, Epr, Dr, Rm, Sm, V)

rm = zeros(1,length(Sm));
[~,L]=qr(Bp-Ep*V*(Epr\(Bpr)), 0);

% Evaluate residual at all sampling points
for ii=1:length(Sm)
    rm(ii)= norm(L*((Dr/Rm)*((Apr-Sm(ii)*Epr)\(Bpr))-eye(size(Bp,2))));
    %     rm2(ii)= norm((Ap-Sm(ii)*Ep)*V*((Apr-Sm(ii)*Epr)\Bpr)-Bp);
end

% Compute maximum
[residual, jz] = max(abs(rm));
shift=Sm(jz);

end

function samples = sampleBorder(EtaSet, Opts)

if Opts.isReal || ~any(imag(EtaSet)) % only sampling points on real axis
    if Opts.log
        samples = logspace(log10(EtaSet(1)), log10(EtaSet(end)), Opts.minSamplePoints);
    else
        samples = linspace(EtaSet(1), EtaSet(end), Opts.minSamplePoints);
    end
else
    samples = EtaSet;
    % Add border points until minBorderpoints is reached
    while length(samples) < Opts.minSamplePoints
        newSamples=(samples(1:end-1) + samples(2:end))/2;
        samples = sort([samples, newSamples], 'ComparisonMethod', 'real');
    end
end

end

function samples = sampleInterior(EtaSet, Opts)

EtaSet = sort(EtaSet, 'ComparisonMethod', 'real');

samples = [];
%Generate mesh
if Opts.log
    x= logspace(log10(min(real(EtaSet)))+1e-12, log10(max(real(EtaSet)))-1e-12, 100);
    y= linspace(0, max(imag(EtaSet))-1e-12, 100);
else
    x= linspace(min(real(EtaSet))+1e-12, max(real(EtaSet))-1e-12, 100);
    y= linspace(0, max(imag(EtaSet))-1e-12, 100);
end

for i=1:length(x)
    % Get neighbors from EtaSet
    set = sort([x(i), EtaSet], 'ComparisonMethod', 'real');
    ind = find(set == x(i));
    if length(ind) == 1
        neighbor = [set(ind-1), set(ind+1)];
        samples = [samples, x(i)+1i*y(y<=imag(neighbor(1))+(imag(neighbor(2))-imag(neighbor(1)))/(real(neighbor(2))-real(neighbor(1)))*(x(i)-real(neighbor(1))))];
    elseif length(ind) == 2 && ~(ind(1)==1 || ind(2) == length(set)) % x is real part of a member of EtaSet
        samples = [samples, x(i)+1i*y(y<=max(imag(set(ind))))];
    else
        % Do not sample border again
    end
end

end

function [rMax, sMin, sMax, lFix, Opts] = parseInputs(sys, rMax, varargin)

% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:trksmPH:wrongInput', 'Original model is not an object of the phs-class');
end

% Catch DAE systems
if sys.isDAE
    error('trksmPH currently only supports pHODE systems.')
end

if rMax >= sys.dim
    error('MORpH:trksmPH:wrongInput', 'Desired reduced order is greater than or equal to original order.');
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.tol = 1e-3;
OptsAdmissible.minSamplePoints = 500;
OptsAdmissible.isReal = {false,true};
OptsAdmissible.border = {true,false};
OptsAdmissible.log = {true,false};
OptsAdmissible.structurePreservation.structurePreservation = 'scaling';
OptsAdmissible.structurePreservation.reorthog = false;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% Remaining parameters
switch length(varargin)
    case 0
        sMin = 1e-5;
        sMax = 1e3;
        lFix = 0;
    case 2
        sMin = varargin{1};
        sMax = varargin{2};
        lFix = 0;
    case 3
        sMin = varargin{1};
        sMax = varargin{2};
        lFix = varargin{3};
    otherwise
        error('MORpH:trksmPH:badInputPattern',...
            'The number of inputs is not correct.')
end

if lFix > size(sys.G,2)
    error('MORpH:trksmPH:wrongInput', 'Number of tangential directions has to be smaller or equal to number of inputs');
end
if real(sMin) >= real(sMax) || real(sMin) <= 0 || real(sMax) <= 0
    error('MORpH:trksmPH:wrongInput', 'Please provide sMin, sMin with positive real parts real(sMin) < real(sMax)');
end

end