function sysr = sobmor(sys, q, varargin)
% SOBMOR - obtaining a H-inf-optimal reduced pH system via parameter optimization
% Syntax:
%   sysr = sobmor(sys, q)
%   sysr = sobmor(sys, q, options)
%
% Description:
%       sysr = sobmor(sys,q, options) returns a reduced PH system with given order q
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:          phs object
%       - q:            desired reduced order
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .maxIter                  Maximum number of iterations 
%                                       (leveled least-squares: max iterates per loop)
%                                       [{1e5} / positive integer]
%           - .adaptiveSampling         Use logarithmic sampling adaption from [2]
%                                       [{true} / false]
%           - .samplePoints             Sample points used for optimization 
%                                       if .adaptiveSampling is false
%                                       [{[logspace(-5,3,800),0,1e-8,1e-7,1e-6,1e4,1e5,1e6]'} / vector]
%           - .additionalSamplePoints	Additional sample points to those 
%                                       from adaptPH used for optimization
%                                       [{logspace(-8,5,14)'} / vector]
%           - .selectionScheme          Selection scheme for gamma
%                                       [{'sequence'} / 'bisection']
%           - .gammaSequence            Gamma sequence values
%                                       [{flip(logspace(-11,2,500)} / vector]
%           - .gammaMax                 Upper bound for gamma in bisection search algorithm
%                                       [{1} / positive double]
%           - .tolBisection             Tolerance for determining if bisection is converged
%                                       [{1e-1} / positive double]
%           - .tolTermination           Tolerance for determining if if leveled least-squares 
%                                       optimization converged to local optimum
%                                       [{1e-14} / positive double]
%           - .tolLsa                   Tolerance for adding new sample point 
%                                       in logarithmic sampling adaption
%                                       [{0.5} / positive double]
%           - .printLevel               Level of detail printed to console regarding optimization progress
%                                       [{0} / 1 / 2 / 3]
%           - .initialReduction         Initial reduction algorithm to compute initial parameter vector
%                                       [{'irkaPH'} / 'adaptPH' / 'cirkaPH']
%           - .plot                     Make intermediate SVD plots including distribution of sample points
%                                       [{false} / true]
%           - .summary                  Print summary at the end
%                                       [{true} / false]
%           - .granso.*                 Other options that will be passed on to the used method;
%                                       Please refer to documentation (doc granso).
%           - .adaptPH.*                Other options that will be passed on to the used method;
%                                       Please refer to documentation (doc adaptPH).
%           - .cirkaPH.*                Other options that will be passed on to the used method;
%                                       Please refer to documentation (doc cirkaPH).
%           - .phs.*                    Other options that will be passed on to the used method;
%                                       Please refer to documentation (doc phs).
%
% Output Arguments:
%       - sysr:     reduced pH system (phsRed object)
%
% See Also:
%       irkaPH, cirkaPH, granso
%
% References:
%       [1] P. Schwerdtner and M. Voigt. Structure Preserving Model Order Reduction 
%           by Parameter Optimization. arXiv Preprint arXiv:2011.07567. 2020. 
%           url: https://arxiv.org/abs/2011.07567.
%       [2] P. Schwerdtner, M. Voigt. Adaptive Sampling for Structure 
%           Preserving Model Order Reduction of Port-Hamiltonian Systems.
%           IFAC-PapersOnLine 54(19) (2021), pp. 143-148.
%       [3] F. E. Curtis, T. Mitchell, and M. L. Overton. “A BFGSSQP method for 
%           nonsmooth, nonconvex, constrained optimization and its evaluation 
%           using relative minimization profiles.” In: Optimization Methods 
%           and Software 32.1 (2017), pp. 148–181.%
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
%               This function uses the software package 'GRANSO' [3], 
%               please refer to its specific license.
%-----------------------------------------------------------------------

%% Start time measurement
startTimeMeasurement = tic;

%% Input parsing
[sys,options] = parseInputs(sys,varargin{:});

% Calculate dimensions for parameter vector
dims.nx = q;
dims.nu = size(sys.G,2);
dims.nJ = dims.nx*(dims.nx-1)/2;
dims.nR = dims.nx*(dims.nx+1)/2;
dims.nQ = dims.nR;
dims.nG = dims.nx*dims.nu;
dims.nvar = dims.nJ + dims.nR + dims.nQ + dims.nG;

%% Initialization
% Run initial reduction method
switch options.initialReduction
    case 'cirkaPH'
        [sysc, V, W, s0, Rt, sysm, s0mTot, nLU] = cirkaPH(sys, q, options.cirkaPH);
        sysc = makeExplicit(sysc);
        samplePoints = imag(s0);
        samplePoints = samplePoints(samplePoints>0);

    case 'irkaPH'
        [sysc, V, s0, Rt, W, nLU] = irkaPH(sys, q, options.irkaPH);
        sysc = makeExplicit(sysc);
        samplePoints = imag(s0);
        samplePoints = samplePoints(samplePoints>0);

    case 'adaptPH'
        [sysc, f, s0] = adaptPH(sys, q, options.adaptPH);
        sysc = makeExplicit(sysc);
        samplePoints = s0;
    case 'rand'
        x0 = rand(dims.nvar,1);
        samplePoints = [];
end

if ~strcmp(options.initialReduction,'none')
    x0 = initializeSystem(sysc);
end

% For GRANSO
options.granso.x0 = x0;

% Initialize sample points
if ~options.adaptiveSampling
    samplePoints = options.samplePoints;
else
    % Bring options.additionalSamplePoints in correct form if necessary
    if isrow(options.additionalSamplePoints)
        options.additionalSamplePoints = options.additionalSamplePoints';
    end

    % Bring samplePoints in correct form if necessary
    if isrow(samplePoints)
        samplePoints = samplePoints';
    end

    % Concatenate sample point vectors
    samplePoints = [samplePoints; options.additionalSamplePoints];
end

% Initialize samples
samples = cell(length(samplePoints),3);
for i=1:length(samplePoints)
    samples{i,1}=samplePoints(i);
end

%% Leveled least-squares
% Evaluate original system at points in samplePoints
for j = 1:length(samplePoints)
    samples{j,2} = evalTF(sys,samples{j,1});
end

% Initialize test sample points
sampleTest = cell(0,2);

% Define figure numbers
np = 10;

switch options.selectionScheme
    case 'sequence'
        gammas = options.gammaSequence;
        % Leveled least-squares iteration
        for j = 1:length(gammas)
            if options.printLevel > 0
                disp(['Gamma: ', num2str(gammas(j))])
            end

            % Check if sampling adaption is used
            if options.adaptiveSampling
                % Get current pH system
                sysr = makephs(options.granso.x0, dims, options);

                % Compute new sampling points
                [samples, sampleTest] = lsa(gammas(j), sys, sysr, samples, sampleTest, options);
            end
            % Declare the objective funciton for GRANSO
            objFun = @(theta)objectiveFunction(theta, samples, gammas(j), dims, options);

            % Start optimization with GRANSO
            solution = granso(dims.nvar, objFun, [], [], options.granso);

            % Initialize x0 of next iteration with optimization result
            options.granso.x0 = solution.best.x;

            if options.plot
                plotdist(sys, sysr, samples, gammas(j), np)
                np=np+1;
            end

            if solution.best.f > options.tolTermination
                final_gamma = gammas(j);
                break
            elseif j == length(gammas)
                final_gamma = gammas(j);
            end
        end

    case 'bisection'
        % Initialize gammaMin and gammaMax
        gammaMax = options.gammaMax;
        gammaMin = 0;
        gamma = (gammaMax+gammaMin)/2;

        while (gammaMax-gammaMin)/(gammaMax+gammaMin)>options.tolBisection
            %while options.tolBisection*gamma<(gammaMax-gammaMin)
            gamma = (gammaMax+gammaMin)/2;
            if options.printLevel > 0
                disp(['Gamma: ', num2str(gamma)])
            end

            % Check if sampling adaption is used
            if options.adaptiveSampling
                % Get current pH system
                sysr = makephs(options.granso.x0, dims, options);

                % Compute new sampling points
                [samples, sampleTest] = lsa(gamma, sys, sysr, samples, sampleTest, options);
            end

            % Declare the objective funciton for GRANSO
            objFun = @(theta)objectiveFunction(theta, samples, gamma, dims, options);

            % Start optimization with GRANSO
            solution = granso(dims.nvar, objFun, [], [], options.granso);

            % Initialize x0 of next iteration with optimization result
            options.granso.x0 = solution.best.x;

            if options.plot
                plotdist(sys, sysr, samples, gamma, np)
                np=np+1;
            end

            % Check whether objective functional is bigger or smaller than current gamma and adjust boundaries accordingly
            if solution.best.f > options.tolTermination
                gammaMin = gamma;
            else
                gammaMax = gamma;
            end
        end
end

% Construct final pH system
sysr = makephs(solution.best.x, dims, options);

% Perform input validation if not done earlier
if ~options.phs.inputValidation
    phs.inputValidation(sysr);
end

% End time measurement for whole procedure
ellapsedtime = toc(startTimeMeasurement);

%% Define output
sysr.parameters = options;
sysr.parameters.samples = samples;
sysr.info.time = ellapsedtime;
sysr.info.dims = dims;
sysr.method = @sobmor;

if strcmp(options.selectionScheme,'bisection')
    sysr.info.gamma = gamma;
    sysr.info.gammaMax = gammaMax;
    sysr.info.gammaMin = gammaMin;
end

if strcmp(options.selectionScheme,'sequence')
    sysr.info.gamma = final_gamma;
    sysr.parameters.gammas = gammas;
end

% Define output of initial reduction method
switch options.initialReduction
    case 'cirkaPH'
        sysr.parameters.cirkaPH.s0 = s0;
        sysr.parameters.cirkaPH.V = V;
        sysr.parameters.cirkaPH.W = W;
        sysr.parameters.cirkaPH.Rt = Rt;
        sysr.parameters.cirkaPH.sysr = sysc;
        sysr.parameters.cirkaPH.sysm = sysm;
        sysr.parameters.cirkaPH.s0mTot = s0mTot;
        sysr.parameters.cirkaPH.nLU = nLU;

    case 'irkaPH'
        sysr.parameters.irkaPH.s0 = s0;
        sysr.parameters.irkaPH.V = V;
        sysr.parameters.irkaPH.W = W;
        sysr.parameters.irkaPH.Rt = Rt;
        sysr.parameters.irkaPH.sysr = sysc;
        sysr.parameters.irkaPH.nLU = nLU;

    case 'adaptPH'
        sysr.parameters.adaptPH.sysr = sysc;
        sysr.parameters.adaptPH.f = f;
        sysr.parameters.adaptPH.s0 = s0;

end

% End sobmor
end

%% ===================== AUXILIARY FUNCTIONS ======================
function [sys,options] = parseInputs(sys,varargin)
narginchk(1,2);

% Check if options is provided
if ~isempty(varargin) && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

if isa(sys,'phs') || isa(sys,'phsRed')
    if sys.isDAE
        error('sobmor does not work with DAE models.')
    end
else
    if ~isempty(sys.E) && svds(sparse(sys.E),1,'smallest') < 1e-12
        error('sobmor currently only supports ODE systems.')
    end
end

if ~(isa(sys,'phs') || isa(sys,'phsRed'))
    options.initialReduction = 'none';
end

% Check if sample points are provided and correct dimensions if necessary
if isfield(options, 'samplePoints')
    if isrow(options.samplePoints)
        options.samplePoints = options.samplePoints';
    end
end

% Check if provided sample points are imaginary
if isfield(options, 'samplePoints')
    if ~isreal(options.samplePoints)
        if sum(real(options.samplePoints)) == 0
            warning('Shift vector might contain real and imaginary parts.')
            options.samplePoints = imag(options.samplePoints);
        else
            error('Shift vector contains real and imaginary parts.')
        end
    end
end

% Admissible value sets
OptsAdmissible.cirkaPH = struct();                                  % Make sure cirkaPH options exists
OptsAdmissible.irkaPH = struct();                                   % Make sure cirkaPH options exists
OptsAdmissible.phs = struct();                                      % Make sure phs options exists
OptsAdmissible.phs.inputValidation = {false,true};                  % Enable input validation
OptsAdmissible.adaptiveSampling = {true,false};                     % Logarithmic sampling adaption from [2]
OptsAdmissible.initialReduction = {'irkaPH','adaptPH','cirkaPH','none'};   % Initial reduction algorithm for initialization
OptsAdmissible.printLevel = {0,1,2,3};                              % Define amount of console output
OptsAdmissible.plot = {false,true};                                 % Make intermediate SVD plots with distribution of sample points
OptsAdmissible.selectionScheme = {'sequence','bisection'};          % Gamma selection scheme

% Default values
OptsAdmissible.cirkaPH.summary = false;                             % No summary at end of cirkaPH
OptsAdmissible.gammaMax = 1;                                        % Maximum cut-off value
OptsAdmissible.tolBisection = 1e-1;                                 % Tolerance for termination of bisection algorithm
OptsAdmissible.tolTermination = 1e-14;                              % Tolerance for termination of optimization
OptsAdmissible.tolLsa = 0.5;                                      	% Tolerance for adding new sample point in logarithmic sampling adaption
OptsAdmissible.maxIter = 1e5;                                       % Maximum number of iterations per optimization step
OptsAdmissible.adaptPH.selectw = 'error';                           % Specify how the new shift is chosen in each iteration -> Hinf error
OptsAdmissible.adaptPH.initw = 'eigs';                              % Specify method for the selection of the initial frequencies -> Eigenvalues
OptsAdmissible.additionalSamplePoints = logspace(-8,5,14)';         % Sepcify sample points additional to those from adaptPH
samplePoints = logspace(-5,3,800);
samplePoints = [samplePoints, 0, 1e-8, 1e-7, 1e-6, 1e4, 1e5, 1e6]'; % Sample points if no adaptive sampling is done
OptsAdmissible.samplePoints = samplePoints;                         % sample points for leveled least-squares
OptsAdmissible.gammaSequence = flip(logspace(-11,2,500));        	% Gamma sequence values
OptsAdmissible.phs.verbose = false;                                 % Neglect console output from phs
OptsAdmissible.cirkaPH.maxIter = 500;                               % Maximum number of outer cirka iterations
OptsAdmissible.cirkaPH.tol = 1e-3;                                  % Convergence tolerance
OptsAdmissible.cirkaPH.irkaPH.maxIter = 500;                        % Maximum number of inner irka iterations
OptsAdmissible.cirkaPH.stopCrit = 's0';                             % Choose shifts as convergence criterium
OptsAdmissible.irkaPH.initShifts = 'zeros';                         % Initialize shifts as zeros
OptsAdmissible.irkaPH.stopCrit = 's0';                              % Choose shifts as convergence criterium
OptsAdmissible.irkaPH.tol = 1e-3;                                   % Convergence tolerance
OptsAdmissible.irkaPH.maxIter = 500;                                % Maximum number of irka iterations
OptsAdmissible.irkaPH.arnoldiPH.structurePreservation = 'specialInverse';   % Specify structure preservation scheme

% Suppress warnings from low-level functions to reduce console output
OptsAdmissible.phs.verbose = false;
OptsAdmissible.cirkaPH.summary = false;
OptsAdmissible.cirkaPH.verbose = false;
OptsAdmissible.irkaPH.summary = false;
OptsAdmissible.irkaPH.verbose = false;

options = phsMOR_parseOpts(options,OptsAdmissible);

% Transfer input values to GRANSO
options.granso.print_level = options.printLevel;
options.granso.maxit = options.maxIter;

thirdPartyCheck('GRANSO');
end

function [cost, grad] = objectiveFunction(theta, samples, gamma, dims, options)
% Get current reduced parameterized pH system
sysr = makephs(theta, dims, options);

% Preallocate singular values and singular vectors
k = size(samples,1);
S = zeros(k, dims.nu);

% Get parts of parameter vector corresponding to R and Q
[~, thetaR, thetaQ, ~] = unziptheta(theta, dims);

% Preallocate gradient
df = zeros(dims.nvar,1);

% Set L to zero
L = 0;

for i = 1:k
    % Get transfer function matrices
    G = samples{i,2};
    Gr = evalTF(sysr, samples{i,1});

    % Calculate singular values
    [u,s,v] = svd(G-full(Gr));
    S(i,:) = diag(s);

    % Calculate F
    F = 1i*samples{i,1}*eye(dims.nx)-(sysr.J-sysr.R)*sysr.Q;

    for j = 1:dims.nu
        if S(i,j)>gamma
            L = L + (S(i,j)-gamma)^2;

            f = S(i,j)-gamma;

            % Get supporting variables
            Y1 = -sysr.Q/F*sysr.G*v(:,j)*u(:,j)'*sysr.G.'*sysr.Q/F;
            Y2 = F\sysr.G*v(:,j)*(u(:,j)'*sysr.G.'+u(:,j)'*sysr.G.'*sysr.Q/F*(sysr.J-sysr.R));

            % Calculate gradients at current sample point
            dJ = real(sutv(-Y1+Y1.'));
            dR = real(utv(vtu(thetaR)*Y1+vtu(thetaR)*Y1.'));
            dQ = real(utv(vtu(thetaQ)*Y2+vtu(thetaQ)*Y2.'));
            dG = real(ftv((v(:,j)*u(:,j)'*sysr.G.'*sysr.Q/F).'+sysr.Q/F*sysr.G*v(:,j)*u(:,j)'));

            % Concatenate and sum gradients
            df = df + f*[dJ; dR; dQ; dG];
        end
    end
end

% Calculate cost
cost = L/gamma;

% Calculate gradient
grad = -df*2/gamma;
end

function [thetaJ, thetaR, thetaQ, thetaG] = unziptheta(theta, dims)
% Decompose parameter vector
thetaJ = theta(1:dims.nJ);
thetaR = theta(dims.nJ+1:dims.nJ+dims.nR);
thetaQ = theta(dims.nJ+dims.nR+1:dims.nJ+dims.nR+dims.nQ);
thetaG = theta(dims.nJ+dims.nR+dims.nQ+1:end);
end

function sys = makephs(theta, dims, options)
[thetaJ, thetaR, thetaQ, thetaG] = unziptheta(theta, dims);

% Construct system matrices
J = vtsu(thetaJ)'-vtsu(thetaJ);
R = vtu(thetaR)'*vtu(thetaR);
Q = vtu(thetaQ)'*vtu(thetaQ);
G = vtf(thetaG, dims.nx, dims.nu);

sys = phsRed(J, R, Q, G, options.phs);
end

function x0 = initializeSystem(sys)
% Decompose system matrices
JS = triu(sys.J,1);
try
    RS = chol(sys.R);
catch
    RS = chol(sys.R+1e-12*eye(size(sys.R)));
end
try
    QS = chol(sys.Q);
catch
    QS = chol(sys.Q+1e-12*eye(size(sys.Q)));
end
GS = sys.G;

% Assign parameter vectors
thJ = sutv(JS);
thR = utv(RS);
thQ = utv(QS);
thG = ftv(GS);
x0 = [thJ;thR;thQ;thG];
end

function [samples, sampleTest] = lsa(gamma, sys, sysr, samples, sampleTest, options)
% Logarithmic sampling adaption for leveled least-squares optimization
n = 1;
k = size(samples,1);

% Preallocate and compute the error at old sample points
for j = 1:k
    samples{j,3} = evalErr(samples{j,1}, samples{j,2}, sysr);
end

while n > 0
    % Sort sample points in ascending order
    samples = sortrows(samples,1);

    n = 0;

    for i = 1:k-1
        % Compute new test sample point
        f = 10^((log10(samples{i,1})+log10(samples{i+1,1}))/2);

        % Check if transfer function at sample point has already been
        % computed
        idxTest = cell2mat(sampleTest(:,1)) == f;

        if all(idxTest==0)
            sampleTest{end+1,2} = evalTF(sys,f);
            sampleTest{end,1} = f;
            G = sampleTest{end,2};
        else
            G = sampleTest{idxTest,2};
        end

        % Compute error at test sample point
        errTest = evalErr(f, G, sysr);

        % Compute difference quotients between old points and new one
        d1 = (errTest-samples{i,3})/(f-samples{i,1});
        d2 = (samples{i+1,3}-errTest)/(samples{i+1,1}-f);

        % Get maximum values
        gammastar = max([samples{i,3},samples{i+1,3}]);
        dstar = max([d1,d2]);

        % Check if test sample point should be included in sample points
        if dstar*(samples{i+1,1}-samples{i,1})>=2*(gammastar+options.tolLsa*gamma)-(samples{i,3}+samples{i+1,3})
            n = n+1;
            samples{end+1,1} = f;
            samples{end,2} = full(G);
            samples{end,3} = errTest;
        end
    end
end
end

function E = evalErr(w, G, sysr)
% Evaluate transfer function at new frequency w
Gr = evalTF(sysr, w);

% Calculate 2-norm of error system
E = norm(G-full(Gr), 2);
end

function G = evalTF(sys, s)
% Evaluate transfer functions at given s
if isa(sys, 'phs')
    G = (sys.G+sys.P)'*sys.Q/(s*1i*sys.E-(sys.J-sys.R)*sys.Q)*(sys.G-sys.P)+(sys.S+sys.N);
else
    G = sys.C/(s*1i*sys.E-sys.A)*sys.B+sys.D;
end
end

function plotdist(sys, sysr, samples, gamma, np)
% Get frequencies
samplePoints = cell2mat(samples{:,1});

% current Hinf norm
[hinf,fpeak] = norm(sys-sysr,Inf);

% Frequency response and singular values at sample points
fs = cell(length(samplePoints),1);
data = zeros(length(samplePoints),1);
for i = 1:length(samplePoints)
    fs{i} = evalTF(sys-sysr, samplePoints(i));
    data(i) = norm(fs{i},2);
end

% Sigmaoptions
plotoptions = sigmaoptions;
plotoptions.MagUnits = 'abs';

% Marker Size
mS = 6;
ml = 'rx--';

% Plot
figure(np)
hold on
sigmaplot(sys-sysr,plotoptions)

for i = 1:length(samplePoints)
    plot([samplePoints(i) samplePoints(i)], [0 data(i)],ml,'MarkerSize',mS);
end

plot([fpeak fpeak], [0 hinf],'gx--','MarkerSize',mS);

yline(gamma,'m--');

end