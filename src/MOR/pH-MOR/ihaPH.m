function [sysr, status, output] = ihaPH(sys, s0, varargin)
% ihaPH - obtaining a H-infinity tuned reduced PH system via feedthrough optimization
%
% Syntax:
%   sysr = ihaPH(sys, redOrder)
%   sysr = ihaPH(sys, s0)
%   sysr = ihaPH(sys, s0, options)
%   sysr = ihaPH(sys, s0, Rt)
%   sysr = ihaPH(sys, s0, Rt, options)
%   [sysr, status, output] = ihaPH(sys,redOrder)
%
% Description:
%       sysr = ihaPH(sys, redOrder) returns a reduced PH system with given order redOrder
%
%       sysr = ihaPH(sys, redOrder, s0, Rt) returns a reduced PH system
%                with given order redOrder, with s0 and Rt as the initial shifts
%                and tangential directions for cirkaPH
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object
%       - redOrder: desired reduced order
%       - s0:       (alternatively) initial guess for shifts for cirkaPH
%       *Optional Input Arguments:*
%       - Rt:       initial guess for right tangential directions for cirkaPH
%       - Opts:  structure with execution parameters
%           - .S0                   initial value for S
%                                   [{rand(m,m)} / symmetric positive definite matrix]
%           - .M0                   initial value for M
%                                   [{rand(m,m)} / skew-symmetric matrix]
%           - .initialPopulation    initial population for gradient free optimization
%                                   [{[]} / cell array containing structs (MIMO) or integers (SISO)]
%           - .optimization         optimization algorithm;
%                                   [{'llsq'} / 'granso' / 'trustregions' / 'bfgsmanopt' /... 
%                                   'arc' / 'steepestdescent' / 'conjugategradient' /... 
%                                   'barzilaiborwein' / 'particleswarm' / 'neldermead']
%           - .hinfnorm.method      method to calculate hinf norm of error system
%                                   [{'matlab'} / 'gsp' / 'sps' / 'sss']
%           - .modelFct             use modelFunction from cirkaPH instead of original system for optimization
%                                   [{0} / 1]
%           - .maxIter              max number of iterations (leveled least-squares: max iterates per loop)
%                                   [{500} / positive integer]
%           - .tolgradnorm          Tolerance for reaching (approximate) optimality/stationarity
%                                   [{1e-8} / positive double]
%           - .adaptiveSampling     use Adaptive Sampling (AS) algorithm [6] for leveled least-squares optimization: 
%                                   0: no AS / 1: AS initialized with user input / ...
%                                   2: AS initialized with adaptPH / 3: AS initialized with (c)irkaPH
%                                   [{3} / 0 / 1 / 2]
%           - .samplePoints       	sample points used for leveled least-squares optimization if .adaptiveSampling is 0 or 1
%                                   [{[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 0 1e1 1e2 1e3 1e4 1e5]'} / vector]
%           - .bisection            use bisection search algorithm from Schwerdtner 2020 [5] for leveled least-squares optimization
%                                   [{false} / true]
%           - .gamma                series of cut-off values for leveled least-squares cost function if .bisection is false
%                                   [{[1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01,..., 0.0000001]} / vector]
%           - .gammaMax             upper bound for gamma in bisection search algorithm
%                                   [{1} / positive double]
%           - .tolBisection         tolerance for determining if bisection is converged
%                                   [{1e-3} / positive double]
%           - .tolTermination       tolerance for determining if if leveled least-squares optimization converged to local optimum
%                                   [{1e-14} / positive double]
%           - .tolLSA               tolerance for adding new sample point in logarithmic sampling adaption
%                                   [{0.5} / positive double]
%           - .finiteDiffOrder      order of magnitude for finite differences approximation of gradients
%                                   [{11} / positive integer]
%           - .printLevel           level of detail printed to console regarding optimization progress
%                                   [{0} / 1 / 2 / 3]
%           - .initialReduction     initial reduction algorithm to compute initial shifts s0
%                                   [{cirkaPH} / irkaPH]
%           - .cmptRealBasis        factorization scheme for computing real basis V and xi
%                                   [{'svd'} / 'qr']
%           - .plot                 make intermediate SVD plots including distribution of sample points
%                                   [{false} / true]
%           - .summary              print summary at the end of ihaPH
%                                   [{true} / false]
%           - .hinfnorm.*           Other options that will be passed on to the third-party 
%                                   functions 'hinorm' or 'linorm_subsp'
%                                   Please refer to documentation (hinorm / linorm_subsp).
%           - .manopt.*             Other options that will be passed on to the third-party toolbox 'Manopt';
%                                   Please refer to documentation (manopt).
%           - .granso.*             Other options that will be passed on to the third-party toolbox 'GRANSO';
%                                   Please refer to documentation (granso).
%           - .adaptPH.*          	Other options that will be passed on to the function 'adaptPH';
%                                   Please refer to documentation (adaptPH).
%           - .cirkaPH.*            Other options that will be passed on to the function 'cirkaPH';
%                                   Please refer to documentation (cirkaPH).
%           - .phs.*                Other options that will be passed on to the class 'phs';
%                                   Please refer to documentation (phs).
%
% Output Arguments:
%       - sysr:     phsRed object, final reduced PH system with optimized feedthrough term
%       - status:   output of optimization procedure
%       - output:   struct containing the output of other functions used
%
% See Also:
%       cirkaPH, irkaPH, granso, linorm_subsp, hinorm, manopt, demo_ihaPH
%
% References:
%       [1] C. Beattie and S. Gugercin. Interpolatory projection methods for
%           structure-preserving model reduction. Systems & Control Letters, 
%           58(3):225–232, 2009.
%       [2] G. Flagg, C. Beattie, and S. Gugercin. Interpolatory H∞
%           model reduction. Systems & Control Letters, 62(7):567–574, 2013.
%       [3] N. Boumal et al. “Manopt, a Matlab Toolbox for Optimization on Manifolds.” 
%           In: Journal of Machine Learning Research 15.42 (2014), pp. 1455–1459. 
%           url: https://www.manopt.org.
%       [4] F. E. Curtis, T. Mitchell, and M. L. Overton. “A BFGSSQP method for 
%           nonsmooth, nonconvex, constrained optimization and its evaluation 
%           using relative minimization profiles.” In: Optimization Methods 
%           and Software 32.1 (2017), pp. 148–181.
%       [5] P. Schwerdtner and M. Voigt. Structure Preserving Model Order Reduction 
%           by Parameter Optimization. arXiv Preprint arXiv:2011.07567. 2020. 
%           url: https://arxiv.org/abs/2011.07567.
%       [6] P. Schwerdtner, M. Voigt. Adaptive Sampling for Structure 
%           Preserving Model Order Reduction of Port-Hamiltonian Systems.
%           IFAC-PapersOnLine 54(19) (2021), pp. 143-148.
%       [7] T. Moser, J. Durmann, and B. Lohmann. “SurrogateBased H2 Model 
%           Reduction of Port-Hamiltonian Systems." In: 2021 Eur. Control 
%           Conf. (ECC). 2021, pp. 2058–2065.
%       [8] P. Benner and M. Voigt. “A structured pseudospectral
%           method for H∞-norm computation of large-scale descriptor
%           systems.” In: Mathematics of Control, Signals, and Systems
%           26.2 (2013), pp. 303–338.
%       [9] P. Schwerdtner and M. Voigt. "Computation of the H∞-Norm Using 
%           Rational Interpolation.” In: IFACPapersOnLine 51.25 (2018). 
%           9th IFAC Symposium on Robust Control Design ROCOND 2018, pp. 84–89.
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
%               If you use any of the software packages 'Manopt' [3], 'GRANSO' [4], 
%               'hinorm' [8] or 'linorm_subsp' [9], please refer to their specific license.
%-----------------------------------------------------------------------

%% Start time measurement for whole procedure
startTimeMeasurement = tic;

%% Input Parsing
[s0, Rt, dims, optFramework, options] = parseInputs(sys, s0, varargin{:});

%% (c)irkaPH
% Start time measurement initial reduction
tRed = tic;

switch options.initialReduction
    case 'cirkaPH'
        [sysc, Vc, Wc, s0, Rt, sysm, s0mTot, nLU] = cirkaPH(sys, s0, Rt, options.cirkaPH);
    case 'irkaPH'
        [sysc, Vc, s0, Rt, Wc, nLU] = irkaPH(sys, dims.redOrder, options.irkaPH);
end

% End time measurement initial reduction
tReduction = toc(tRed);

%% IHA
% Compute real basis and xi
k = 1;
Vreell = zeros(sys.dim,length(s0));
Rtreell = zeros(size(sys.G,2),length(s0));

% Make V and Rt real
while(k < length(s0)+1)
    v = ((s0(k)*sys.E)-(sys.J-sys.R)*sys.Q)\(sys.G-sys.P)*Rt(:,k);
    if(abs(imag(s0(k))) > 0)
        Vreell(:,k) = real(v);
        Rtreell(:,k) = real(Rt(:,k));
        Vreell(:,k+1) = imag(v);
        Rtreell(:,k+1) = imag(Rt(:,k));
        k = k + 2;
    else
        Vreell(:,k) = real(v);
        Rtreell(:,k) = real(Rt(:,k));
        k = k + 1;
    end
end

% Choose factorization method
switch options.cmptRealBasis
    case 'svd'
        [V,Do,Uo] = svd(Vreell,0);  % Vreell = V*Do*Uo' --> Tw = Do*Uo'
        RtRe = Rtreell/(Do*Uo');

    case 'qr'
        [V,Tv2] = qr(Vreell,0);
        RtRe = Rtreell/Tv2;
end

% Compute new system sysn with V and W
Qn = V'*sys.Q*V;
W = sys.Q*V/Qn;
Jn = W'*sys.J*W;
Rn = W'*sys.R*W;
Gn = W'*sys.G;
Pn = W'*sys.P;
if isdiag(sys.E) && all(abs(diag(sys.E) - 1) < 1e-12)
    En = eye(size(Jn));
else
    En = W'*sys.E*V;
end
sysn = phs(Jn, Rn, Qn, Gn, En, Pn, sysc.S, sysc.N, options.phs);

% Check if sysn is pH-system
if ~options.phs.inputValidation
    phs.inputValidation(sysn);
end

% Compute xi
xi = sysn.Q\RtRe';

% Define whether model function or original system is used for optimization
switch options.modelFct
    case 0
        system = sys;
    case 1
        system = sysm;
end

% Make original system sss
syssss = phs2sss(system);

% Initialize number of iterations for optimization
iter = 0;

% Apply correct optimization framework
switch optFramework
    % Granso
    case 'granso'
        % Specify problem for GRANSO
        objFun = @(theta)objFunDirect(theta, syssss, sysn, xi, dims, options);

        % Start time measurement optimization
        toptim = tic;

        solution = granso(dims.nvar,objFun,[],[],options.granso);

        % End time measurement optimization
        tOptimization = toc(toptim);

        % Adjust output from granso
        [Sopt, Nopt] = getOptimParams(solution.best.x, dims, options);

        % Fetch output
        hinf = solution.best.f;
        status = solution;
        iter = solution.iters;

        % Manopt
    case 'manopt'
        % Specify problem for manopt
        if dims.m > 1
            if strcmp(options.optimization, 'neldermead') || strcmp(options.optimization, 'particleswarm')
                problem.M = productmanifold(struct('S', euclideanfactory(dims.ns,1), ...
                    'M', euclideanfactory(dims.nm,1)));
            else
                problem.M = productmanifold(struct('S', sympositivedefinitefactory(dims.m), ...
                    'M', skewsymmetricfactory(dims.m)));
            end
        else
            if strcmp(options.optimization, 'neldermead') || strcmp(options.optimization, 'particleswarm')
                problem.M = euclideanfactory(dims.m);
            else
                problem.M = positivefactory(dims.m);
            end
        end

        % Define cost function
        problem.cost = @(theta)costFunDirect(theta, syssss, sysn, xi, dims, options);

        % Disable warning because of missing gradient and Hessian
        if ~strcmp(options.optimization, 'neldermead') && ~strcmp(options.optimization, 'particleswarm')
            warning('off', 'manopt:getGradient:approx')
            warning('off', 'manopt:getHessian:approx')
        end

        % Start time measurement optimization
        toptim = tic;

        % Select optimization algorithm
        switch options.optimization
            case 'trustregions'
                [solution, hinf, status, ~] = trustregions(problem, options.manopt.x0, options.manopt);
            case 'arc'
                [solution, hinf, status, ~] = arc(problem, options.manopt.x0, options.manopt);
            case 'steepestdescent'
                [solution, hinf, status, ~] = steepestdescent(problem, options.manopt.x0, options.manopt);
            case 'conjugategradient'
                [solution, hinf, status, ~] = conjugategradient(problem, options.manopt.x0, options.manopt);
            case 'barzilaiborwein'
                [solution, hinf, status, ~] = barzilaiborwein(problem, options.manopt.x0, options.manopt);
            case 'bfgsmanopt'
                [solution, hinf, status, ~] = rlbfgs(problem, options.manopt.x0, options.manopt);
            case 'particleswarm'
                [solution, hinf, status, ~] = pso(problem, options.initialPopulation, options.manopt);
            case 'neldermead'
                [solution, hinf, status, ~] = neldermead(problem, options.initialPopulation, options.manopt);
        end

        % End time measurement optimization
        tOptimization = toc(toptim);

        % Adjust output from manopt
        [Sopt, Nopt] = getOptimParams(solution, dims, options);
        iter = status(end).iter;

        % Leveled least-squares
    case 'llsq'
        % Get sample points
        switch options.adaptiveSampling
            case 2
                % Start time measurement adaptPH
                tadaptPH = tic;

                % AdaptPH to get initial sample points
                [sysa, hinfAdaptPH, samplePointsAdaptPH] = adaptPH(system,dims.redOrder,options.adaptPH);

                % End time measurement adaptPH
                tAdaptPH = toc(tadaptPH);

                % Add additional, initial sample points
                samplePoints = [samplePointsAdaptPH,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5]';

            case 3
                if isrow(s0)
                    samplePoints = imag(s0');
                else
                    samplePoints = imag(s0);
                end
                samplePoints = samplePoints(samplePoints>0);
                samplePoints = [samplePoints;1e-8;1e-7;1e-6;1e-5;1e-4;1e-3;1e-2;1e-1;1e0;1e1;1e2;1e3;1e4;1e5];

            otherwise
                samplePoints = options.samplePoints;
                if isrow(samplePoints)
                    samplePoints = samplePoints';
                end

                % End switch-case Adaptive Sampling
        end

        % Number of sample points before adaptiveSampling
        if options.adaptiveSampling~=0
            numInitialSamples = length(samplePoints);
        end

        % Initialize samples
        samples = cell(length(samplePoints),3);
        for i=1:length(samplePoints)
            samples{i,1}=samplePoints(i);
        end

        % Evaluate sample points and measure time
        tevalSamples = tic;
        samples = evaluateSamplePoints(system,samples);
        tEvaluation = toc(tevalSamples);

        % Initialize test sample point cell array and array
        sampleTest = cell(0,2);

        % Start time measurement optimization
        toptim = tic;

        % Check if bisection algorithm or gamma sequence is used
        if ~options.bisection
            % Get gamma sequence values
            gamma = options.gamma;

            % Leveled least-squares iteration
            for j = 1:length(gamma)
                [solution, samples, sampleTest, iter, options] = runLLSQ(system, sysn, xi, gamma(j), samples, sampleTest, iter, dims, options);
                if solution.best.f > options.tolTermination
                    g = gamma(j);
                    break
                end
            end

            % Check if optimization ran out of gamma values
            ise = evalin( 'base', 'exist(''g'',''var'') ~= 0' );
            if ise
                g = inf;
            end

        else % Gamma bisection
            gammaMax = options.gammaMax;
            gammaMin = 0;
            gammas = cell(1);
            counter = 1;
            np = 10;
            while (gammaMax-gammaMin)/(gammaMax+gammaMin)>options.tolBisection
                gamma = (gammaMax+gammaMin)/2;
                gammas{counter} = gamma;
                [solution, samples, sampleTest, iter, options] = runLLSQ(system, sysn, xi, gamma, samples, sampleTest, iter, dims, options);

                if options.plot
                    plotdist(sys, sysn, xi, gamma, solution, samplePoints, np, dims, options)
                    np = np+1;
                end

                if solution.best.f > options.tolTermination
                    gammaMin = gamma;
                else
                    gammaMax = gamma;
                end
                counter = counter+1;
            end

            % End if Bisection
        end

        % End time measurement optimization
        tOptimization = toc(toptim);

        % Adjust output from granso
        [Sopt, Nopt] = getOptimParams(solution.best.x, dims, options);

        % Fetch output
        status = solution;

        % End switch-case optimization framework
end

% Define pH System with feedthrough matrix (S+M)
[Jr, Rr, Qr, Gr, Er, Pr] = chooseSystemMatrices(sysn, xi, Sopt, Nopt);

% Make reduced system
sysr = phsRed(Jr, Rr, Qr, Gr, Er, Pr, sysn.S+Sopt, sysn.N+Nopt, options.phs);

% Check if sysr is port-Hamiltonian, if not already checked
if ~options.phs.inputValidation
    phs.inputValidation(sysr);
end

% End time measurement for whole procedure
ellapsedtime = toc(startTimeMeasurement);

%% Set output
% Calculate hinf norm for leveled least-squares or when modelFct is
% used with original system
if options.modelFct == 1
    syssss = phs2sss(sys);
    sysrsss = phs2sss(sysr);
    if options.printNorm
        hinf = calculateNorm(syssss, sysrsss, options);
    end
elseif strcmp(optFramework, 'llsq')
    sysrsss = phs2sss(sysr);
    if options.printNorm
        hinf = calculateNorm(syssss, sysrsss, options);
    end
end

% Calculate approximate hinf norm when modelFct is used
if options.modelFct == 1
    sysmsss = phs2sss(sysm);
    if options.printNorm
        hinfModelFct = calculateNorm(sysmsss, sysrsss, options);
    end
end

% Calculate norm of (c)irkaPH
sssCirkaPH = phs2sss(sysc);
if options.printNorm
    hinfCirkaPH = calculateNorm(syssss, sssCirkaPH, options);
end

% Display a warning if (c)irkaPH solves the problem better
if options.printNorm
    if hinfCirkaPH < hinf
        warning('ROM from (c)irkaPH achieves smaller Hinf norm.')
    end

    % Improvement in percent over (c)irkaPH
    perImprovement = (1-(hinf/hinfCirkaPH))*100;
end

% Display summary
if options.summary
    fprintf("======================================\n");
    fprintf(strcat("IHAPH stopped after ",num2str(iter)," iterations\n"));
    if options.printNorm
        fprintf(strcat("Hinf norm of error system: ",num2str(hinf),"\n"));
        fprintf(strcat("Improvement over (C)IRKAPH: ",num2str(round(perImprovement,2)),"%%\n"));
    end
    fprintf(strcat("Elapsed time: ",num2str(round(ellapsedtime,2)),"s\n"));
    fprintf("======================================\n");
end

% Define output for sysr
if options.printNorm
    sysr.info.hinf = hinf;
    sysr.info.improvement = perImprovement;
end
sysr.info.dims = dims;
sysr.info.tOptimization = tOptimization;
sysr.info.tReduction = tReduction;
sysr.info.time = ellapsedtime;
sysr.method = @ihaPH;
sysr.parameters = options;
sysr.parameters.xi = xi;
sysr.parameters.sysn = sysn;
if options.printNorm
    if options.modelFct == 1
        sysr.info.hinfModelFct = hinfModelFct;
    end
end

switch options.initialReduction
    case 'cirkaPH'
        % Define output for cirkaPH
        output.cirkaPH.sysr = sysc;
        output.cirkaPH.V = Vc;
        output.cirkaPH.W = Wc;
        output.cirkaPH.s0 = s0;
        output.cirkaPH.Rt = Rt;
        output.cirkaPH.sysm = sysm;
        output.cirkaPH.s0mTot = s0mTot;
        output.cirkaPH.nLU = nLU;
        if options.printNorm
            output.cirkaPH.hinf = hinfCirkaPH;
        end

    case 'irkaPH'
        % Define output for irkaPH
        output.irkaPH.sysr = sysc;
        output.irkaPH.V = Vc;
        output.irkaPH.W = Wc;
        output.irkaPH.s0 = s0;
        output.irkaPH.Rt = Rt;
        output.irkaPH.nLU = nLU;
        if options.printNorm
            output.irkaPH.hinf = hinfCirkaPH;
        end
end

% Output number of optimization steps
if strcmp(options.optimization,'llsq') || strcmp(options.optimization,'granso')
    status.iterations = iter;
end

% Define output from adaptPH if adaptiveSampling is done
if options.adaptiveSampling == 2 && strcmp(options.optimization, 'llsq')
    output.adaptPH.sysa = sysa;
    output.adaptPH.hinf = hinfAdaptPH;
    output.adaptPH.samplePoints = samplePointsAdaptPH;
    output.tAdaptPH = tAdaptPH;
end

% Define output for leveled least-squares
if strcmp(optFramework, 'llsq')
    sysr.parameters.samplePoints = samplePoints;
    sysr.info.tEvaluation = tEvaluation;

    % Compute number of added samples from adaptive sampling
    if options.adaptiveSampling
        sysr.info.numAdaptiveSamples = length(samplePoints)-numInitialSamples;
    end

    if options.bisection
        sysr.info.gamma = gamma;
        sysr.info.gammaMax = gammaMax;
        sysr.info.gammaMin = gammaMin;
        sysr.parameters.gammas = gammas;
    else
        sysr.info.gamma = g;
    end
end

end

%% ===================== AUXILIARY FUNCTIONS ======================
function [s0, Rt, dims, optFramework, options] = parseInputs(sys, s0, varargin)
narginchk(2,4);

% Check phs & phsRed input type
if ~(isa(sys,'phs') || isa(sys,'phsRed'))
    error('MORpH:ihaPH:wrongInput', 'Original model is not an object of the phs-class.');
end

if sys.isDAE
    error('ihaPh currently only supports pHODE systems.')
end

% Check if input s0 is reduced order (redOrder) or vector of frequencies
if length(s0) == 1 && mod(s0,1) == 0 && s0 ~= 0
    % Reduced order specified
    dims.redOrder = s0;
    s0 = zeros(1,s0);
else
    dims.redOrder = length(s0);
end

% Check if options is provided
if ~isempty(varargin) && isstruct(varargin{end})
    options = varargin{end};
    varargin = varargin(1:end-1);
else
    options = struct();
end

% Check remaining inputs
switch length(varargin)
    case 0
        Rt = ones(size(sys.G,2),length(s0));
    case 1
        Rt = varargin{1};
end

% Original system dimension
dims.n = size(sys.Q,1);

% Input Dimension
dims.m = size(sys.G,2);

% Parameter numbers
dims.ns = dims.m*(dims.m+1)/2;
dims.nm = dims.m*(dims.m-1)/2;
dims.nvar = dims.ns + dims.nm;

% Set hinfnorm calculation method if not specified
if ~isfield(options,'hinfnorm')
    if dims.n < 500
        options.hinfnorm.method = 'matlab';
    else
        options.hinfnorm.method = 'gsp';
    end
else
    if ~isfield(options.hinfnorm,'method')
        if dims.n < 500
            options.hinfnorm.method = 'matlab';
        else
            options.hinfnorm.method = 'gsp';
        end
    end
end

% Admissible value sets
OptsAdmissible.optimization = {'llsq','granso','trustregions','bfgsmanopt','arc','steepestdescent','conjugategradient','barzilaiborwein','particleswarm','neldermead'}; % Optimization algorithms
OptsAdmissible.hinfnorm.method = {'matlab','gsp', 'sps','sss'}; % Hinf norm computation methods
OptsAdmissible.printNorm = {false,true};                % Define if the Hinf norms of fom and rom are calculated and compared
OptsAdmissible.printLevel = {0,1,2,3};                  % Level of detail printed to console regarding optimization progress
OptsAdmissible.cirkaPH = struct();                      % Make sure cirkaPH options exists
OptsAdmissible.irkaPH = struct();                       % Make sure irkaPH options exists
OptsAdmissible.phs = struct();                          % Make sure phs options exists
OptsAdmissible.granso = struct();                       % Make sure granso options exists
OptsAdmissible.manopt = struct();                       % Make sure manopt options exists
OptsAdmissible.adaptPH = struct();                      % Make sure adaptPH options exists
OptsAdmissible.modelFct = {0,1};                        % Use model function from cirkaPH for optimization
OptsAdmissible.phs.inputValidation = {false,true};      % Validate if system is port-Hamiltonian
OptsAdmissible.adaptiveSampling = {3,0,1,2};            % Use adaptive sampling
OptsAdmissible.bisection = {false,true};                % Use gamma bisection search
OptsAdmissible.summary = {true,false};                  % Print summary at end of ihaPH
OptsAdmissible.initialReduction = {'cirkaPH','irkaPH'}; % Initial reduction algorithm to compute initial shifts s0
OptsAdmissible.plot = {false,true};                     % Make intermediate SVD plots with distribution of sample points
OptsAdmissible.cmptRealBasis = {'svd','qr'};            % Factorization scheme for computation of real basis V and xi

% Default values
OptsAdmissible.maxIter = 500;                           % Maximum number of optimization iterations
OptsAdmissible.finiteDiffOrder = 11;                    % Order of finite differences scheme used for optimization with granso
OptsAdmissible.tolgradnorm = 1e-8;                      % Tolerance for reaching (approximate) optimality/stationarity
OptsAdmissible.tolLSA = 0.5;                        	% Tolerance for adding new sample point in logarithmic sampling adaption
OptsAdmissible.adaptPH.selectw = 'error';               % Choose next shift as the maximizer of H-inf norm of the error model in adaptPH
OptsAdmissible.adaptPH.initw = 'eigs';                  % Initial frequencies are first chosen in a wide range and then reinitialized in adaptPH
OptsAdmissible.cirkaPH.maxIter = 500;                   % Maximum number of outer cirka iterations
OptsAdmissible.cirkaPH.tol = 1e-6;                      % Convergence tolerance
OptsAdmissible.cirkaPH.irkaPH.maxIter = 500;            % Maximum number of inner irka iterations
OptsAdmissible.cirkaPH.stopCrit = 's0';                 % Choose shifts as convergence criterium
OptsAdmissible.irkaPH.initShifts = 'zeros';             % Initialize shifts as zeros
OptsAdmissible.irkaPH.stopCrit = 's0';                  % Choose shifts as convergence criterium
OptsAdmissible.irkaPH.tol = 1e-6;                       % Convergence tolerance
OptsAdmissible.irkaPH.maxIter = 500;                    % Maximum number of irka iterations
OptsAdmissible.irkaPH.arnoldiPH.structurePreservation = 'specialInverse';   % Specify structure preservation scheme

% Default values that are checked in second routine
OptsAdmissible.samplePoints = [];                     	% samplePoints for leveled least-squares cost function
OptsAdmissible.gamma = [];                              % Cut-off values for leveled least-squares cost function

% Suppress warnings from low-level functions to reduce console output
OptsAdmissible.phs.verbose = false;
OptsAdmissible.cirkaPH.summary = false;
OptsAdmissible.cirkaPH.verbose = false;
OptsAdmissible.irkaPH.summary = false;
OptsAdmissible.irkaPH.verbose = false;

% Check if initialPopulation is provided
if isfield(options, 'initialPopulation')
    if isa(options.initialPopulation, 'cell')
        if sys.isMIMO
            if numel(options.initialPopulation) ~= dims.nvar+1
                error('Initial population: wrong dimension')
            else
                for k = 1:numel(options.initialPopulation)
                    vertex = options.initialPopulation{k};
                    if numel(vertex.S) ~= dims.ns || numel(vertex.M) ~= dims.nm
                        error('Initial population: wrong dimension')
                    end
                end
            end
        else
            if numel(options.initialPopulation) ~= 2
                error('Initial population: wrong dimension')
            end
        end
    else
        error('Initial population: wrong data type')
    end
else
    % No initial population -> provide manopt empty array
    options.initialPopulation = [];
end

options = phsMOR_parseOpts(options,OptsAdmissible);

% If irkaPH is used as initial reduction method, the modelFct can't be selected
if strcmp(options.initialReduction, 'irkaPH')
    options.modelFct = 0;
end

% Set default values for hinfnorm computation methods
switch options.hinfnorm.method
    case 'matlab'
        OptsAdmissible.hinfnorm.tol = 0.01;                             % Relative accuracy of the H∞ norm
    case 'gsp'
        OptsAdmissible.hinfnorm.boydbalak = {1,0};                      % Use matlabs norm command in Schwerdtner 2018 Hinf norm computation (0->AB13HD only on linux)
        OptsAdmissible.hinfnorm.doLoewner = {0};                        % Specify that the Loewner framework is not used
        OptsAdmissible.hinfnorm.fnHandle = {0};                         % Specify that no function handle is passed
        OptsAdmissible.hinfnorm.initialPoints = linspace( 0, 100, 10);  % Initial frequencies for Hinf norm computation in Schwerdtner 2018
    otherwise
        OptsAdmissible.hinfnorm = struct();                             % Make sure options.hinfnorm exists
end

% Set flag according to optimization framework
if strcmp(options.optimization, 'llsq')
    flag_opt = 'granso';
    optFramework = 'llsq';
elseif strcmp(options.optimization, 'granso')
    flag_opt = 'granso';
    optFramework = 'granso';
else
    flag_opt = 'manopt';
    optFramework = 'manopt';
end
options.flag_opt = flag_opt;

switch flag_opt
    % Set granso specific options
    case 'granso'
        OptsAdmissible.granso.print_level = options.printLevel;     % Level of detail printed to console regarding optimization progress
        OptsAdmissible.granso.maxit = options.maxIter;              % Max number of iterations
        OptsAdmissible.granso.opt_tol = options.tolgradnorm;        % Tolerance for reaching (approximate) optimality/stationarity

        % Check options only necessary for leveled least-squares cost function
        if strcmp(options.optimization, 'llsq')
            % Check if sample points are provided if no adaptive sampling is used
            if isempty(options.samplePoints) && options.adaptiveSampling ~= 2 && options.adaptiveSampling ~= 3
                options.samplePoints = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 0 1e1 1e2 1e3 1e4 1e5]';
            end

            % Check if gammas are provided if no bisection search is used
            if isempty(options.gamma) && ~options.bisection
                options.gamma = [1, 0.5, 0.2, 0.1, 0.05,0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002,0.0001, 0.00005, 0.00002, 0.00001, 0.000005,0.000002, 0.000001, 0.0000005, 0.0000002,0.0000001];
            end
            OptsAdmissible.tolBisection = 1e-3;                     % Bisection tolerance for leveled least-squares
            OptsAdmissible.tolTermination = 1e-14;                	% Termination tolerance for leveled least-squares
            OptsAdmissible.gammaMax = 1;                        	% Maximal gamma value for bisection search
        end

        % Check if options.granso.x0 is provided
        if isfield(options.granso, 'x0')
            if ~sys.isMIMO
                options.S0 = options.granso.x0;
            else
                options.S0 = options.granso.x0.S;
                options.M0 = options.granso.x0.M;
            end
        end

        % Define initial values for S and M
        % Check if S0 is provided
        if isfield(options,'S0')
            % Check if S0 is provided as parameter vector with dims.ns entries
            if isvector(options.S0) && numel(options.S0) == dims.ns
                S0 = options.S0;
                % Check dimension of S0 as matrix
            elseif size(options.S0,1) == dims.m && size(options.S0,2) == dims.m
                % Check positive-definiteness of S0
                try
                    % Get initial parameter vector for S0
                    S0 = utv(chol(options.S0));
                catch ME
                    error('S0 must be symmetric and positive definite.')
                end
            else
                error('Matrix dimensions wrong: S0')
            end
        else
            % Random Initialization if S0 is not provided
            S0 = 0.0001*randn(dims.ns,1);
            % Make sure S0 is positive if SISO
            if ~sys.isMIMO
                S0 = abs(S0);
            end
        end

        if sys.isMIMO
            % Check if M0 is provided
            if isfield(options,'M0')
                % Check dimension of M0
                % Check if M0 is provided as parameter vector with dims.nm entries
                if isvector(options.M0) && numel(options.M0) == dims.nm
                    M0 = options.M0;
                elseif size(options.M0,1) == dims.m && size(options.M0,2) == dims.m
                    % Check skew-symmetry of M0
                    if ~issymmetric(options.M0,'skew')
                        error('M0 must be skew-symmetric.')
                    end
                    % Make entries of M0 a vector
                    M0 = sutv(options.M0);
                else
                    error('Matrix dimensions wrong: M0')
                end
            else
                % Random initialization if M0 is not provided
                M0 = zeros(dims.nm,1);
            end
        else
            M0 = [];
        end

        % Check if column or row vector
        if isrow(S0)
            S0 = S0';
        end
        if isrow(M0)
            M0 = M0';
        end
        options.granso.x0 = [S0;M0];

        % Set manopt specific options
    case 'manopt'
        OptsAdmissible.manopt.verbosity = options.printLevel;           % Level of detail printed to console regarding optimization progress
        OptsAdmissible.manopt.maxiter = options.maxIter;                % Max number of iterations
        OptsAdmissible.manopt.tolgradnorm = options.tolgradnorm;        % Tolerance for reaching (approximate) optimality/stationarity

        % Check if options.manopt.x0 is provided
        if isfield(options.manopt, 'x0') && ~strcmp(options.optimization, 'neldermead') && ~strcmp(options.optimization, 'particleswarm')
            if ~sys.isMIMO
                options.S0 = options.manopt.x0;
            else
                options.S0 = options.manopt.x0.S;
                options.M0 = options.manopt.x0.M;
            end
        end

        if ~strcmp(options.optimization, 'neldermead') && ~strcmp(options.optimization, 'particleswarm')
            if isfield(options,'S0')
                % Check dimension of S0
                if size(options.S0,1) == dims.m && size(options.S0,2) == dims.m
                    % Check positive-definiteness of S0
                    try
                        chol(options.S0);
                    catch ME
                        error('S0 must be symmetric and positive definite.')
                    end
                    % Check if SISO or MIMO
                    if sys.isMIMO
                        options.manopt.x0.S = options.S0;
                    else
                        options.manopt.x0 = options.S0;
                    end
                    % Check if S0 is provided as parameter vector with dims.ns entries
                elseif isvector(options.S0) && numel(options.S0) == dims.ns
                    S0 = vtu(options.S0)'*vtu(options.S0);
                    % Check if SISO or MIMO
                    if sys.isMIMO
                        options.manopt.x0.S = S0;
                    else
                        options.manopt.x0 = S0;
                    end
                else
                    error('Matrix dimensions wrong: S0')
                end
            else
                % Random Initialization if S0 is not provided
                if sys.isMIMO
                    K = randn(dims.ns,1);
                    H = vtu(K);
                    options.manopt.x0.S = 0.0001*(H*H');
                else
                    options.manopt.x0 = abs(0.0001*randn(1));
                    options.M = [];
                end
            end

            if sys.isMIMO
                % Check if M0 is provided
                if isfield(options,'M0')
                    % Check dimension of M0
                    if size(options.M0,1) == dims.m && size(options.M0,2) == dims.m
                        % Check skew-symmetry of M0
                        if ~issymmetric(options.M0,'skew')
                            error('M0 must be skew-symmetric.')
                        end
                        options.manopt.x0.M = options.M0;
                        % Check if M0 is provided as parameter vector with dims.nm entries
                    elseif isvector(options.M0) && numel(options.M0) == dims.nm
                        options.manopt.x0.M = vtsu(options.M0)-vtsu(options.M0)';
                    else
                        error('Matrix dimensions wrong: M0')
                    end
                else
                    options.manopt.x0.M = zeros(dims.m);
                end
            end
        end

        % End switch-case Adjust optimization options
end

options = phsMOR_parseOpts(options,OptsAdmissible);

switch optFramework
    case 'manopt'
        thirdPartyCheck('Manopt');
    otherwise
        thirdPartyCheck('GRANSO');
end
if options.printNorm
    switch options.hinfnorm.method
        case 'gsp'
            thirdPartyCheck('LINORM');
        case 'sps'
            thirdPartyCheck('HINORM');
    end
end
end

function C = costFunDirect(theta, sys, sysc, xi, dims, options)
% Cost Function for direct optimization
% Objective Function for direct optimization with GRANSO
sysr = makephs(theta, dims, sysc, xi, options);

% Make reduced system sss
sysr = phs2sss(sysr);

% Calculate the H-inf norm
C = calculateNorm(sys, sysr, options);
end

function [cost, grad] = objFunDirect(theta, sys, sysc, xi, dims, options)
% Objective Function for direct optimization with GRANSO
% Function handle for finite differences scheme
fun = @(theta)costFunDirect(theta, sys, sysc, xi, dims, options);

% Cost function evaluation
cost = costFunDirect(theta, sys, sysc, xi, dims, options);

% Gradient evaluation
grad = finiteDiff(theta, dims, options.finiteDiffOrder, fun);
end

function [sol, samples, sampleTest, iter, options] = runLLSQ(system, sysn, xi, gamma, samples, sampleTest, iter, dims, options)
% Check if adaptive sampling is active
if options.adaptiveSampling
    % Get pH system with current parametrization for adaptive Sampling
    sysr = makephs(options.granso.x0, dims, sysn, xi, options);

    % Perform adaptive sampling
    [samples, sampleTest] = lsa(gamma, system, sysr, samples, sampleTest, options);
end

% Declare objective function for GRANSO
objFun = @(theta)objFunLLSQ(theta, samples, gamma, dims, sysn, xi, options);

% Start optimization with GRANSO
sol = granso(dims.nvar,objFun,[],[],options.granso);

% Make best guess new initial starting point
options.granso.x0 = sol.best.x;
iter = iter + sol.iters;
end

function [cost, grad] = objFunLLSQ(theta, samples, gamma, dims, sysn, xi, options)
% Objective Function for leveled least-squares optimization with GRANSO
% Get pH system with current parametrization
sysr = makephs(theta, dims, sysn, xi, options);

% Preallocate singular values and singular vectors
k = size(samples,1);
S = zeros(k, dims.m);

% Get parameter vector for S
[thetaS, ~] = unzipTheta(theta, dims);

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
    F = 1i*samples{i,1}*sysr.E-(sysr.J-sysr.R)*sysr.Q;

    for j = 1:dims.m
        if S(i,j)>gamma
            L = L + (S(i,j)-gamma)^2;

            f = S(i,j)-gamma;

            % Compute all summands of the gradient
            M1 = xi.'*sysr.Q/(F)*(sysr.G-sysr.P)*v(:,j)*u(:,j)';
            M2 = xi.'*sysr.Q/(F)*(sysr.G-sysr.P)*v(:,j)*u(:,j)'*(sysr.G+sysr.P).'*sysr.Q/(F)*xi;
            M3 = v(:,j)*u(:,j)'*(sysr.G+sysr.P).'*sysr.Q/(F)*xi;
            M4 = v(:,j)*u(:,j)';

            % Compute gradients
            dM = real(sutv(-M1+M1.')+sutv(-M2+M2.')-sutv(-M3+M3.')-sutv(-M4+M4.'));
            dS = real(utv(vtu(thetaS)*(M1+M2-M3-M4)+vtu(thetaS)*(M1+M2-M3-M4).'));

            % Concatenate and sum gradients
            df = df + f*[dS;dM];
        end
    end
end

% Calculate cost
cost = L/gamma;

% Calculate gradient
grad = df*2/gamma;
end

function G = evalTF(sys, s)
% Evaluate transfer functions at given s
if isa(sys, 'phs')
    %G = (sys.G+sys.P)'*sys.Q/(s*1i*eye(size(sys.Q))-(sys.J-sys.R)*sys.Q)*(sys.G-sys.P)+(sys.S+sys.N);
    G = (sys.G+sys.P)'*sys.Q/(s*1i*sys.E-(sys.J-sys.R)*sys.Q)*(sys.G-sys.P)+(sys.S+sys.N);
else
    %G = sys.C/(s*1i*eye(size(sys.A))-sys.A)*sys.B;
    G = sys.C/(s*1i*sys.E-sys.A)*sys.B + sys.D;
end
end

function samples = evaluateSamplePoints(system,samples)
% Evaluate original system at points in samplePoints
for j = 1:size(samples,1)
    samples{j,2} = evalTF(system,samples{j,1});
end
end

function [gradf] = finiteDiff(x0, dims, k, fun)
% Finite differences
% Preallocation of the gradient
gradf = zeros(dims.nvar,1);

% Get function value at x0
f0 = fun(x0);

% Preallocate vector with deflected function values
f1 = zeros(dims.nvar,1);

% Calculate difference quotients for all directions
for i=1:dims.nvar
    % Save initial value
    x1 = x0;
    z = zeros(dims.nvar,1);

    % Declare stepsize
    z(i) = 10^(-k);

    % Make small step in current direction
    x1 = x1 + z;

    % Evaluate function at new point
    f1(i) = fun(x1);

    % Compute gradient with differential quotient
    gradf(i) = (f1(i)-f0)/10^(-k);
end
end

function hinf = calculateNorm(syssss, sysrsss, options)
switch options.hinfnorm.method
    % sss toolbox method
    case 'sss'
        % Call sss toolbox Hinf norm calculation routine
        hinf = norm(syssss-sysrsss, Inf, options.hinfnorm);

        % Matlab built-in
    case 'matlab'
        % Convert sss objects to ss objects
        sys = ss(syssss);
        sysr = ss(sysrsss);

        % Call matlabs Hinf norm calculation routine
        hinf = norm(sys-sysr, Inf, options.hinfnorm.tol);

        % Greedy subspace method
    case 'gsp'
        % Get error system
        errorsys = getErrorsystem(syssss, sysrsss);

        % Call Hinf norm calculation routine from Schwerdtner 2018
        hinf = linorm_subsp(errorsys, options.hinfnorm);

        % Structured pseudospectral method
    case 'sps'
        % Get error system
        errorsys = getErrorsystem(syssss, sysrsss);

        % Call Hinf norm calculation routine from Benner 2014
        hinf = hinorm(errorsys.E, errorsys.A, errorsys.B, errorsys.C, full(errorsys.D), options.hinfnorm);
end
end

function [S, M] = getOptimParams(theta, dims, options)
% Granso
if strcmp(options.flag_opt, 'granso')
    % Get parameter vectors for S and M
    [thetaS, thetaM] = unzipTheta(theta, dims);

    % Form S from upper triangular matrix
    GS = vtu(thetaS);
    S = GS'*GS;

    % Form M from strictly upper matrix
    GM = vtsu(thetaM);
    M = GM-GM';

    % Gradient-free optimization
elseif strcmp(options.optimization, 'neldermead') || strcmp(options.optimization, 'particleswarm')
    if dims.m > 1
        % Form M from strictly upper matrix
        GS = vtu(theta.S);
        GM = vtsu(theta.M);
        M = GM-GM';
    else
        % For SISO M = 0 and theta contains only S
        GS = vtu(theta);
        M = 0;
    end
    % Form S from upper triangular matrix
    S = GS'*GS;

    % Manopt with gradient
else
    if dims.m > 1
        % Get S and M from theta
        S = theta.S;
        M = theta.M;
    else
        % SISO: thata contains only S
        S = theta;
        M = 0;
    end
end
end

function [thetaS, thetaN] = unzipTheta(theta, dims)
thetaS = theta(1:dims.ns);
thetaN = theta(dims.ns+1:end);
end

function errorsys = getErrorsystem(sys, sysr)
% Make error system of type ss
errorsystem = sys-sysr;

% Convert error system matrices to sparse matrices - not of type ss
errorsys.A = errorsystem.A;
errorsys.B = errorsystem.B;
errorsys.C = errorsystem.C;
errorsys.D = errorsystem.D;
errorsys.E = errorsystem.E;
end

function sysr = makephs(theta, dims, sys, xi, options)
% Get S and N from theta
[S, N] = getOptimParams(theta, dims, options);

% Choose the system matrices with new feedtrough terms S and N
[J, R, Q, G, E, P] = chooseSystemMatrices(sys, xi, S, N);

% Get port-Hamiltonian system
sysr = phs(J, R, Q, G, E, P, sys.S+S, sys.N+N, options.phs);
end

function [Js, Rs, Qs, Gs, Es, Ps] = chooseSystemMatrices(sys, xi, S, M)
% Compute the system matrices with the new feedtrough terms S and M
Js = sys.J - xi*M*xi.';
Rs = sys.R + xi*S*xi.';
Qs = sys.Q;
Gs = sys.G+xi*M;
Es = sys.E;
Ps = sys.P-xi*S;
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
        if dstar*(samples{i+1,1}-samples{i,1})>=2*(gammastar+options.tolLSA*gamma)-(samples{i,3}+samples{i+1,3})
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

if issparse(Gr)
    G = full(G);
    Gr = full(Gr);
end

% Calculate 2-norm of error system
E = norm(G-Gr, 2);
end

function plotdist(sys, sysn, xi, gamma, solution, samplePoints, np, dims, options)
% Adjust output from granso
[Sopt, Nopt] = getOptimParams(solution.best.x, dims, options);

% Define pH System with feedthrough matrix (S+M)
[Jr, Rr, Qr, Gr, Er, Pr] = chooseSystemMatrices(sysn, xi, Sopt, Nopt);

% Make reduced system
sysr = phsRed(Jr, Rr, Qr, Gr, Er, Pr, sysn.S+Sopt, sysn.N+Nopt, options.phs);

% current Hinf norm
[hinf,fpeak] = norm(sys-sysr,Inf);

% Frequency response and singular values at sample points
fs = cell(length(samplePoints),1);
data = zeros(length(samplePoints),1);
for i = 1:length(samplePoints)
    fs{i} = freqresp(sys-sysr, samplePoints(i));
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