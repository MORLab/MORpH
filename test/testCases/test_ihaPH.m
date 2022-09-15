classdef test_ihaPH < matlab.unittest.TestCase
    % TEST test_ihaPH - Tests the functionalities of ihaPH
    %
    % Press <F5> or enter "runtests("test_ihaPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'ihaPH' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Output: Check if each additional output can be obtained
    %       - Options: Check if additional options work correctly
    %       - Functionality: Try more complex reduction to simulate real-world
    %                        application
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

    properties (TestParameter)
        sysSISOn40 = {setup_MassSpringDamperSystem(40,4,4,1)};   % standard MSD system for most of the tests
        sysSISOn300 = {setup_MassSpringDamperSystem(300,4,4,1)};   % order 300
        sysMIMOn40 = {setup_MassSpringDamperSystem(40,4,4,1,'MIMO')}; % MIMO MSD system for MIMO testing
        sysMIMOn300 = {setup_MassSpringDamperSystem(300,4,4,1,'MIMO')}; % order 300
    end

    properties
        tolerance = 1e-7;
    end

    %% Setup / teardown of test environment

    methods (TestClassSetup) % once before all tests
        function lessOutput(testCase)
            % disable warnings to show cleaner output
            warning('off','all');
            testCase.addTeardown(@warning,'on','all');
        end
        function includeSystemSetupFunctions(testCase)
            % make sure that system folder with system setups is included
            if ~(exist('setup_MassSpringDamperSystem.m','file') == 2)
                error("You need to include ""...\MORpH\demos\test_systems"" in your MATLAB path to run the tests!")
            end
        end
    end

    methods (TestClassTeardown) % once after all tests

    end

    methods (TestMethodSetup) % before each test

    end

    methods (TestMethodTeardown) % after each test

    end

    %% TESTS

    methods (Test, TestTags = {'Input'})
        % Test input combinations

        function standardInput(test,sysSISOn40)
            % Very basic input: system and reduced order only
            redSys = ihaPH(sysSISOn40,4);
            test.verifyEqual(redSys.dim, 4);
            test.verifyEqual(redSys.method, @ihaPH);
        end

        %% Optimizer
        % SISO
        function OptimizerLLSQGranso(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = 1e-3;
            opts.optimization = 'llsq';
            [redSys, status, ~] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'llsq');
            testCase.verifyTrue(isfield(status, 'quadprog_failure_rate'));
            testCase.verifyTrue(isfield(redSys.info, 'gamma'));
        end

        function OptimizerDirectGranso(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = 1e-3;
            opts.optimization = 'granso';
            [redSys, status, ~] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'granso');
            testCase.verifyTrue(isfield(status, 'quadprog_failure_rate'));
            testCase.verifyTrue(~isfield(redSys.info, 'gamma'));
            testCase.verifyEqual(redSys.parameters.granso.x0, opts.S0);
        end

        function OptimizerManoptGradient(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = 1e-3;
            opts.optimization = 'trustregions';
            [redSys, status, ~] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'trustregions');
            testCase.verifyTrue(isfield(status, 'Delta'));
            testCase.verifyEqual(redSys.parameters.manopt.x0, opts.S0);
        end

        function OptimizerManoptGradientFree(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.initialPopulation = {0,2};
            opts.optimization = 'neldermead';
            [redSys, status ,~] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'neldermead');
            testCase.verifyTrue(~isfield(status, 'Delta'));
            testCase.verifyTrue(isfield(status, 'costevals'));
        end

        % MIMO - S and M as vector
        function OptimizerDirectGransoMIMO(testCase,sysMIMOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = [1e-3,1e-3,1e-3];
            opts.M0 = 1e-5;
            opts.optimization = 'granso';
            [redSys, status, ~] = ihaPH(sysMIMOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'granso');
            testCase.verifyTrue(isfield(status, 'quadprog_failure_rate'));
            testCase.verifyTrue(~isfield(redSys.info, 'gamma'));
            testCase.verifyEqual(redSys.parameters.granso.x0, [opts.S0';opts.M0]);
        end

        function OptimizerManoptGradientMIMO(testCase,sysMIMOn40)
            % Pass opts struct
            opts.summary = false;
            S = [1e-06,1e-06;1e-06,2e-06];
            M = [0,0.0001;-0.0001,0];
            opts.S0 = [1e-3,1e-3,1e-3];
            opts.M0 = 1e-4;
            opts.optimization = 'trustregions';
            [redSys, status, ~] = ihaPH(sysMIMOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'trustregions');
            testCase.verifyTrue(isfield(status, 'Delta'));
            testCase.verifyEqual(redSys.parameters.manopt.x0.S, S);
            testCase.verifyEqual(redSys.parameters.manopt.x0.M, M);
        end

        function OptimizerManoptGradientFreeMIMO(testCase,sysMIMOn40)
            % Pass opts struct
            opts.summary = false;
            S1.M = 1e-3;
            S2.M = 1e-10;
            S3.M = 1e-2;
            S4.M = 0;
            S5.M = 1e-3;
            % S Must be column vector!
            S1.S = [1;1;0.1];
            S2.S = [0;1;0];
            S3.S = [.5;0;0];
            S4.S = [0;0;1];
            S5.S = [0;.5;1];
            opts.initialPopulation = {S1,S2,S3,S4,S5};
            opts.optimization = 'neldermead';
            [redSys, status, ~] = ihaPH(sysMIMOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'neldermead');
            testCase.verifyTrue(~isfield(status, 'Delta'));
            testCase.verifyTrue(isfield(status, 'costevals'));
        end

        % MIMO - S and M as matrices
        function OptimizerDirectGransoMIMOInitMatrices(testCase,sysMIMOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = [1e-06,1e-06;1e-06,2e-06];
            opts.M0 = [0,0.0001;-0.0001,0];
            S = [1e-3,1e-3,1e-3];
            M = 1e-4;
            opts.optimization = 'granso';
            [redSys, status, ~] = ihaPH(sysMIMOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'granso');
            testCase.verifyTrue(isfield(status, 'quadprog_failure_rate'));
            testCase.verifyTrue(~isfield(redSys.info, 'gamma'));
            testCase.verifyEqual(redSys.parameters.granso.x0, [S';M]);
        end

        function OptimizerManoptGradientMIMOInitMatrices(testCase,sysMIMOn40)
            % Pass opts struct
            opts.summary = false;
            opts.S0 = [1e-06,1e-06;1e-06,2e-06];
            opts.M0 = [0,0.0001;-0.0001,0];
            opts.optimization = 'trustregions';
            [redSys, status, ~] = ihaPH(sysMIMOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.optimization, 'trustregions');
            testCase.verifyTrue(isfield(status, 'Delta'));
            testCase.verifyEqual(redSys.parameters.manopt.x0.S, opts.S0);
            testCase.verifyEqual(redSys.parameters.manopt.x0.M, opts.M0);
        end

        %% Hinf norm methods
        function HinfNormMatlab(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.hinfnorm.method = 'matlab';
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.hinfnorm.method, opts.hinfnorm.method);
        end

        function HinfNormGreedySubspace(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.hinfnorm.method = 'gsp';
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.hinfnorm.method, opts.hinfnorm.method);
        end

        function HinfNormSPS(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.hinfnorm.method = 'sps';
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.hinfnorm.method, opts.hinfnorm.method);
        end

        function HinfNormSSS(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.hinfnorm.method = 'sss';
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.hinfnorm.method, opts.hinfnorm.method);
        end

        %% Gamma update schemes
        function Sequence(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.optimization = 'llsq';
            opts.bisection = false;
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyGreaterThan(length(redSys.parameters.gamma),1);
        end

        function Bisection(testCase,sysSISOn40)
            % Pass opts struct
            opts.summary = false;
            opts.optimization = 'llsq';
            opts.bisection = true;
            redSys = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEmpty(redSys.parameters.gamma);
            testCase.verifyTrue(isfield(redSys.info,'gammaMax'));
        end

        %% StandardInputWithOpts
        function standardInputWithOpts(testCase,sysSISOn40)
            % Pass opts struct
            opts.initialReduction = 'irkaPH';
            opts.summary = false;
            opts.tolBisection = 1e-2;
            opts.tolLSA = 0.7;
            opts.tolTermination = 1e-12;
            opts.printLevel = 0;
            opts.summary = false;
            opts.maxIter = 1000;
            opts.gammaMax = 100;
            opts.bisection = true;
            opts.cmptRealBasis = 'qr';
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            [redSys, ~, out] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.tolBisection,opts.tolBisection)
            testCase.verifyEqual(redSys.parameters.tolLSA,opts.tolLSA)
            testCase.verifyEqual(redSys.parameters.tolTermination,opts.tolTermination)
            testCase.verifyEqual(redSys.parameters.printLevel,opts.printLevel)
            testCase.verifyEqual(redSys.parameters.summary,opts.summary)
            testCase.verifyEqual(redSys.parameters.maxIter,opts.maxIter)
            testCase.verifyLessThan(redSys.info.gammaMax,opts.gammaMax)
            testCase.verifyEqual(redSys.parameters.cmptRealBasis,opts.cmptRealBasis)
            testCase.verifyEqual(out.irkaPH.sysr.parameters.arnoldiPH.structurePreservation,...
                opts.irkaPH.arnoldiPH.structurePreservation)
        end

        function MIMO(testCase,sysMIMOn40)
            % Test MIMO system reduction (with opts)
            opts.initialReduction = 'irkaPH';
            opts.summary = false;
            opts.tolBisection = 1e-2;
            opts.tolLSA = 0.7;
            opts.tolTermination = 1e-12;
            opts.printLevel = 0;
            opts.summary = false;
            opts.maxIter = 1000;
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            [redSys, ~, out] = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.parameters.tolBisection,opts.tolBisection)
            testCase.verifyEqual(redSys.parameters.tolLSA,opts.tolLSA)
            testCase.verifyEqual(redSys.parameters.tolTermination,opts.tolTermination)
            testCase.verifyEqual(redSys.parameters.printLevel,opts.printLevel)
            testCase.verifyEqual(redSys.parameters.summary,opts.summary)
            testCase.verifyEqual(redSys.parameters.maxIter,opts.maxIter)
            testCase.verifyEqual(out.irkaPH.sysr.parameters.arnoldiPH.structurePreservation,...
                opts.irkaPH.arnoldiPH.structurePreservation)
        end

    end

    methods (Test, TestTags = {'Output'})
        % Check if optional output parameters are provided correctly

        function standardOutputGammaSequence(testCase,sysSISOn40)
            % Check if the output object redSys provides additional information
            % in struct parameters (see phsRed class)
            opts.printNorm = true;
            opts.summary = false;
            [redSys, ~,out] = ihaPH(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.method, @ihaPH);
            testCase.verifyEqual(length(out.cirkaPH.s0), 4);

            % property "parameters" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.parameters,'optimization'));
            testCase.verifyTrue(isfield(redSys.parameters,'hinfnorm'));
            testCase.verifyTrue(isfield(redSys.parameters,'initialReduction'));
            testCase.verifyTrue(isfield(redSys.parameters,'cmptRealBasis'));
            testCase.verifyTrue(isfield(redSys.parameters,'modelFct'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptiveSampling'));
            testCase.verifyTrue(isfield(redSys.parameters,'plot'));
            testCase.verifyTrue(isfield(redSys.parameters,'printLevel'));
            testCase.verifyTrue(isfield(redSys.parameters,'bisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolTermination'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolBisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolLSA'));
            testCase.verifyTrue(isfield(redSys.parameters,'maxIter'));
            testCase.verifyTrue(isfield(redSys.parameters,'irkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'phs'));
            testCase.verifyTrue(isfield(redSys.parameters,'cirkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'initialPopulation'));
            testCase.verifyTrue(isfield(redSys.parameters,'granso'));
            testCase.verifyTrue(isfield(redSys.parameters,'manopt'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'summary'));
            testCase.verifyTrue(isfield(redSys.parameters,'finiteDiffOrder'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolgradnorm'));
            testCase.verifyTrue(isfield(redSys.parameters,'samplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'gamma'));
            testCase.verifyTrue(isfield(redSys.parameters,'flag_opt'));
            testCase.verifyTrue(isfield(redSys.parameters,'xi'));
            testCase.verifyTrue(isfield(redSys.parameters,'sysn'));

            % property "info" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.info,'hinf'));
            testCase.verifyTrue(isfield(redSys.info,'dims'));
            testCase.verifyTrue(isfield(redSys.info,'gamma'));
            testCase.verifyTrue(isfield(redSys.info,'improvement'));
            testCase.verifyTrue(isfield(redSys.info,'time'));
        end

        function standardOutputBisection(testCase,sysSISOn40)
            % Check if the output object redSys provides additional information
            % in struct parameters (see phsRed class)
            opts.printNorm = true;
            opts.summary = false;
            opts.bisection = true;
            [redSys, ~, out] = ihaPH(sysSISOn40,4, opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.method, @ihaPH);
            testCase.verifyEqual(length(out.cirkaPH.s0), 4);

            % property "parameters" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.parameters,'optimization'));
            testCase.verifyTrue(isfield(redSys.parameters,'hinfnorm'));
            testCase.verifyTrue(isfield(redSys.parameters,'initialReduction'));
            testCase.verifyTrue(isfield(redSys.parameters,'cmptRealBasis'));
            testCase.verifyTrue(isfield(redSys.parameters,'modelFct'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptiveSampling'));
            testCase.verifyTrue(isfield(redSys.parameters,'plot'));
            testCase.verifyTrue(isfield(redSys.parameters,'printLevel'));
            testCase.verifyTrue(isfield(redSys.parameters,'bisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolTermination'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolBisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolLSA'));
            testCase.verifyTrue(isfield(redSys.parameters,'maxIter'));
            testCase.verifyTrue(isfield(redSys.parameters,'irkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'phs'));
            testCase.verifyTrue(isfield(redSys.parameters,'cirkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'initialPopulation'));
            testCase.verifyTrue(isfield(redSys.parameters,'granso'));
            testCase.verifyTrue(isfield(redSys.parameters,'manopt'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'summary'));
            testCase.verifyTrue(isfield(redSys.parameters,'finiteDiffOrder'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolgradnorm'));
            testCase.verifyTrue(isfield(redSys.parameters,'samplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'gamma'));
            testCase.verifyTrue(isfield(redSys.parameters,'flag_opt'));
            testCase.verifyTrue(isfield(redSys.parameters,'xi'));
            testCase.verifyTrue(isfield(redSys.parameters,'sysn'));

            % property "info" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.info,'hinf'));
            testCase.verifyTrue(isfield(redSys.info,'dims'));
            testCase.verifyTrue(isfield(redSys.info,'gamma'));
            testCase.verifyTrue(isfield(redSys.info,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.info,'gammaMin'));
            testCase.verifyTrue(isfield(redSys.info,'improvement'));
            testCase.verifyTrue(isfield(redSys.info,'time'));
        end

    end

    methods (Test, TestTags = {'Options'})
        % Test optionalal computation arguments (opts struct)
        function MIMOWithBisection(testCase,sysMIMOn40)
            % Turn on bisection
            opts.summary = false;
            opts.gammaMax = 100;
            opts.bisection = true;
            redSys = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyLessThan(redSys.info.gammaMax,opts.gammaMax)
            testCase.verifyEqual(redSys.parameters.bisection,true)
        end

        function MIMOWithGammaSequence(testCase,sysMIMOn40)
            % Turn on bisection
            opts.summary = false;
            opts.gammaMax = 100;
            opts.bisection = false;
            redSys = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyTrue(~isfield(redSys.info,'gammaMax'))
            testCase.verifyEqual(redSys.parameters.bisection,false)
        end

        function MIMOWithIrkaPH(testCase,sysMIMOn40)
            % Use irkaPH as initial reduction method
            opts.summary = false;
            opts.initialReduction = 'irkaPH';
            [redSys, ~,out] = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyTrue(isfield(out,'irkaPH'))
            testCase.verifyEqual(redSys.parameters.initialReduction, opts.initialReduction)
        end

        function MIMOAdaptiveSamplingIrkaPH(testCase,sysMIMOn40)
            % Use predefined samplePoints
            opts.summary = false;
            opts.initialReduction = 'irkaPH';
            opts.adaptiveSampling = 3;
            [redSys,~,out] = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);

            if isrow(out.irkaPH.s0)
                sp = imag(out.irkaPH.s0');
            else
                sp = imag(out.irkaPH.s0);
            end
            sp = sp(sp>0);
            for i = 1:length(sp)
                testCase.verifyTrue(ismember(sp(i),redSys.parameters.samplePoints))
            end
        end

        function MIMOAdaptiveSamplingAdaptPH(testCase,sysMIMOn40)
            % Use predefined samplePoints
            opts.summary = false;
            opts.adaptiveSampling = 2;
            [redSys,~,out] = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);

            for i = 1:length(out.adaptPH.samplePoints)
                testCase.verifyTrue(ismember(out.adaptPH.samplePoints(i),redSys.parameters.samplePoints))
            end
        end

        function MIMOAdaptiveSamplingUserInput(testCase,sysMIMOn40)
            % Use predefined samplePoints
            opts.summary = false;
            opts.adaptiveSampling = 1;
            opts.samplePoints = [1e-3,1e-2,1e-1,1,10];
            redSys = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            for i = 1:length(opts.samplePoints)
                testCase.verifyTrue(ismember(opts.samplePoints(i),redSys.parameters.samplePoints))
            end
        end

        function MIMOWithoutAdaptiveSampling(testCase,sysMIMOn40)
            % Turn off adaptive sampling
            opts.summary = false;
            opts.adaptiveSampling = 0;
            opts.samplePoints = [0.1,1,10,100];
            redSys = ihaPH(sysMIMOn40,2,opts);
            testCase.verifyEqual(redSys.dim, 2);
            testCase.verifyEqual(length(redSys.parameters.samplePoints), 4)
        end

        function MIMOGransoOptions(testCase,sysMIMOn40)
            % Use fixed sample points
            opts.summary = false;
            opts.granso.step_tol = 1e-10;
            opts.granso.maxit = 100;
            redSys = ihaPH(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);

            % ihaPH overwrites max iterations input to granso
            testCase.verifyEqual(redSys.parameters.maxIter, 500);
            testCase.verifyEqual(redSys.parameters.granso.step_tol, opts.granso.step_tol);
        end

    end

    methods (Test, TestTags = {'Functionality'})
        % Simulate more realistic scenarios with larger original models
        % Compare results with other functions (irkaPH)
        function accuracy(testCase, sysSISOn300)
            % make sure reduced system approximates the original
            opts.summary = false;
            opts.printNorm = true;
            [redSys, ~, out] = ihaPH(sysSISOn300,10,opts);
            testCase.verifyGreaterThan(redSys.info.improvement,0);

            s = out.cirkaPH.s0;
            r = out.cirkaPH.Rt;
            diff = zeros(size(s));
            for i = 1:length(s)
                resp_sys = (sysSISOn300.G+sysSISOn300.P)'*sysSISOn300.Q*((s(i)*eye(size(sysSISOn300.Q))-(sysSISOn300.J-sysSISOn300.R)*sysSISOn300.Q)\(sysSISOn300.G-sysSISOn300.P))*r(:,i)+(sysSISOn300.S+sysSISOn300.N)*r(:,i);
                resp_sysr = (redSys.G+redSys.P)'*redSys.Q*((s(i)*eye(size(redSys.Q))-(redSys.J-redSys.R)*redSys.Q)\(redSys.G-redSys.P))*r(:,i)+(redSys.S+redSys.N)*r(:,i);
                diff(i) = norm(resp_sys - resp_sysr);
            end

            maxdiff = max(diff);
            testCase.verifyLessThan(maxdiff,1e-12);
        end

        %% IrkaPH comparison
        function ihaPHVsIrkaPH(testCase, sysSISOn300)
            % Compare results of ihaPH to those of irkaPH
            optsirka.summary = false;
            optsiha.summary = false;
            optsiha.initialReduction = 'irkaPH';
            redSys_irka = irkaPH(sysSISOn300,10,optsirka);
            redSys_iha = ihaPH(sysSISOn300,10,optsiha);
            err_irka = norm(redSys_irka-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            err_iha = norm(redSys_iha-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            testCase.verifyLessThan(err_iha, err_irka);
        end

        function ihaPHVsIrkaPHMIMO(testCase, sysMIMOn300)
            % Compare results of ihaPH to those of irkaPH (MIMO case)
            optsirka.summary = false;
            optsiha.S0 = [1e-3;1e-3;1e-3];
            optsiha.M0 = 0.0001;
            optsiha.summary = false;
            optsiha.initialReduction = 'irkaPH';
            redSys_irka = irkaPH(sysMIMOn300,10,optsirka);
            redSys_iha = ihaPH(sysMIMOn300,10,optsiha);
            err_irka = norm(redSys_irka-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            err_iha = norm(redSys_iha-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            testCase.verifyLessThanOrEqual(err_iha, err_irka);
        end

        %% CirkaPH comparison
        function ihaPHVsCirkaPH(testCase, sysSISOn300)
            % Compare results of ihaPH to those of irkaPH
            optscirka.summary = false;
            optsiha.summary = false;
            optsiha.S0 = 1e-5;
            optsiha.initialReduction = 'cirkaPH';
            redSys_cirka = irkaPH(sysSISOn300,10,optscirka);
            redSys_iha = ihaPH(sysSISOn300,10,optsiha);
            err_cirka = norm(redSys_cirka-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            err_iha = norm(redSys_iha-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            testCase.verifyLessThan(err_iha, err_cirka);
        end

        function ihaPHVsCirkaPHMIMO(testCase, sysMIMOn300)
            % Compare results of ihaPH to those of irkaPH (MIMO case)
            optscirka.summary = false;
            optsiha.S0 = [1e-3;1e-3;1e-3];
            optsiha.M0 = 0.0001;
            optsiha.summary = false;
            optsiha.initialReduction = 'cirkaPH';
            redSys_cirka = irkaPH(sysMIMOn300,10,optscirka);
            redSys_iha = ihaPH(sysMIMOn300,10,optsiha);
            err_cirka = norm(redSys_cirka-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            err_iha = norm(redSys_iha-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            testCase.verifyLessThanOrEqual(err_iha, err_cirka);
        end

        %% Effort Constraint Method comparison
        function ihaPHVsECM(testCase, sysMIMOn300)
            % Compare results of ihaPH to those of Effort Constraint Method (MIMO case)
            redSys_ecm = ecm(sysMIMOn300,10);
            opts.summary = false;
            redSys_iha = ihaPH(sysMIMOn300,10,opts);
            err_ecm = norm(redSys_ecm-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            err_iha = norm(redSys_iha-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            testCase.verifyLessThanOrEqual(err_iha, err_ecm);
        end
    end
end