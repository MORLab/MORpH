classdef test_sobmor < matlab.unittest.TestCase
    % TEST test_sobmor - Tests the functionalities of sobmor
    %
    % Press <F5> or enter "runtests("test_sobmor")" to run this testscript
    %
    % Description:
    %   This script tests the function 'sobmor' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Output: Check if each additional output can be obtained
    %       - Options: Check if additional options work correctly
    %       - Functionality: Try more complex reduction to simulate real-world
    %                       application
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
        sysSISOn40 = {setup_MassSpringDamperSystem(40,2,1,1)};   % standard MSD system for most of the tests
        sysSISOn300 = {setup_MassSpringDamperSystem(300,2,1,1)};   % order 300
        sysMIMOn40 = {setup_MassSpringDamperSystem(40,2,1,1,'MIMO')}; % MIMO MSD system for MIMO testing
        sysMIMOn300 = {setup_MassSpringDamperSystem(300,2,1,1,'MIMO')}; % order 300
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
            redSys = sobmor(sysSISOn40,4);
            test.verifyEqual(redSys.dim, 4);
            test.verifyEqual(redSys.method, @sobmor);
        end

        function standardInputWithOpts(testCase,sysSISOn40)
            % Pass opts struct
            opts.tolBisection = 1e-2;
            opts.tolLsa = 0.7;
            opts.tolTermination = 1e-12;
            opts.printLevel = 0;
            opts.summary = false;
            opts.maxIter = 10;
            opts.gammaMax = 100;
            opts.selectionScheme = 'bisection';
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            redSys = sobmor(sysSISOn40,4,opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.parameters.tolBisection,opts.tolBisection)
            testCase.verifyEqual(redSys.parameters.tolLsa,opts.tolLsa)
            testCase.verifyEqual(redSys.parameters.tolTermination,opts.tolTermination)
            testCase.verifyEqual(redSys.parameters.printLevel,opts.printLevel)
            testCase.verifyEqual(redSys.parameters.summary,opts.summary)
            testCase.verifyEqual(redSys.parameters.maxIter,opts.maxIter)
            testCase.verifyLessThan(redSys.info.gammaMax,opts.gammaMax)
            testCase.verifyEqual(redSys.parameters.selectionScheme,opts.selectionScheme)
            testCase.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation,...
                opts.irkaPH.arnoldiPH.structurePreservation)
        end

        function MIMO(testCase,sysMIMOn40)
            % Test MIMO system reduction (with opts)
            opts.tolBisection = 1e-2;
            opts.tolLsa = 0.7;
            opts.tolTermination = 1e-12;
            opts.printLevel = 0;
            opts.summary = false;
            opts.maxIter = 10;
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.parameters.tolBisection,opts.tolBisection)
            testCase.verifyEqual(redSys.parameters.tolLsa,opts.tolLsa)
            testCase.verifyEqual(redSys.parameters.tolTermination,opts.tolTermination)
            testCase.verifyEqual(redSys.parameters.printLevel,opts.printLevel)
            testCase.verifyEqual(redSys.parameters.summary,opts.summary)
            testCase.verifyEqual(redSys.parameters.maxIter,opts.maxIter)
            testCase.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation,...
                opts.irkaPH.arnoldiPH.structurePreservation)
        end

    end

    methods (Test, TestTags = {'Output'})
        % Check if optional output parameters are provided correctly

        function standardOutputGammaSequence(testCase,sysSISOn40)
        % Check if the output object redSys provides additional information
        % in struct parameters (see phsRed class)
            redSys = sobmor(sysSISOn40,4);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.method, @sobmor);
            testCase.verifyEqual(length(redSys.parameters.irkaPH.s0), 4);

            % property "parameters" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.parameters,'initialReduction'));
            testCase.verifyTrue(isfield(redSys.parameters,'maxIter'));
            testCase.verifyTrue(isfield(redSys.parameters,'selectionScheme'));
            testCase.verifyTrue(isfield(redSys.parameters,'irkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'cirkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'phs'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptiveSampling'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaSequence'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammas'));
            testCase.verifyTrue(isfield(redSys.parameters,'printLevel'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolTermination'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolBisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolLsa'));
            testCase.verifyTrue(isfield(redSys.parameters,'additionalSamplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'samplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'granso'));

            % property "info" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.info,'time'));
            testCase.verifyTrue(isfield(redSys.info,'dims'));
            testCase.verifyTrue(isfield(redSys.info,'gamma'));
        end

        function standardOutputBisection(testCase,sysSISOn40)
            % Check if the output object redSys provides additional information
            % in struct parameters (see phsRed class)
            opts.selectionScheme = 'bisection';
            opts.maxIter = 10;
            redSys = sobmor(sysSISOn40,4, opts);
            testCase.verifyEqual(redSys.dim, 4);
            testCase.verifyEqual(redSys.method, @sobmor);
            testCase.verifyEqual(length(redSys.parameters.irkaPH.s0), 4);

            % property "parameters" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.parameters,'initialReduction'));
            testCase.verifyTrue(isfield(redSys.parameters,'maxIter'));
            testCase.verifyTrue(isfield(redSys.parameters,'selectionScheme'));
            testCase.verifyTrue(isfield(redSys.parameters,'irkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'cirkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'phs'));
            testCase.verifyTrue(isfield(redSys.parameters,'adaptiveSampling'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaSequence'));
            testCase.verifyTrue(isfield(redSys.parameters,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.parameters,'printLevel'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolTermination'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolBisection'));
            testCase.verifyTrue(isfield(redSys.parameters,'tolLsa'));
            testCase.verifyTrue(isfield(redSys.parameters,'additionalSamplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'samplePoints'));
            testCase.verifyTrue(isfield(redSys.parameters,'granso'));

            % property "info" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.info,'time'));
            testCase.verifyTrue(isfield(redSys.info,'dims'));
            testCase.verifyTrue(isfield(redSys.info,'gamma'));
            testCase.verifyTrue(isfield(redSys.info,'gammaMax'));
            testCase.verifyTrue(isfield(redSys.info,'gammaMin'));
        end

    end

    methods (Test, TestTags = {'Options'})
        % Test optionalal computation arguments (opts struct)

        function MIMOWithoutAdaptiveSampling(testCase,sysMIMOn40)
        % Turn off adaptive sampling
            opts.maxIter = 10;
            opts.adaptiveSampling = false;
            opts.samplePoints = [0.1,1,10,100];
            redSys = sobmor(sysMIMOn40,2,opts);
            testCase.verifyEqual(redSys.dim, 2);
            testCase.verifyEqual(length(redSys.parameters.samplePoints), 4)
        end

        function MIMOWithBisection(testCase,sysMIMOn40)
        % Turn on bisection 
            opts.maxIter = 10;
            opts.gammaMax = 100;
            opts.selectionScheme = 'bisection';
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyLessThan(redSys.info.gammaMax,opts.gammaMax)
            testCase.verifyEqual(redSys.parameters.selectionScheme,opts.selectionScheme)
        end

        function MIMOWithCirkaPH(testCase,sysMIMOn40)
        % Use cirkaPH as initial reduction method
            opts.maxIter = 10;
            opts.initialReduction = 'cirkaPH';
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.parameters.cirkaPH.sysr.method,@cirkaPH)
            testCase.verifyEqual(redSys.parameters.initialReduction, opts.initialReduction)
        end

        function MIMOWithAdaptPH(testCase,sysMIMOn40)
        % Use adaptPH as initial reduction method
            opts.maxIter = 10;
            opts.initialReduction = 'adaptPH';
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.parameters.adaptPH.sysr.method,@adaptPH)
            testCase.verifyEqual(redSys.parameters.initialReduction, opts.initialReduction)
        end

        function MIMOWithAdditionalSamplePoints(testCase,sysMIMOn40)
        % Use predefined samplePoints
            opts.maxIter = 10;
            opts.additionalSamplePoints = [12,24,48];
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.method, @sobmor);

            samples = cell2mat(redSys.parameters.samples(:,1));
            testCase.verifyEqual(ismember(12,samples),true);
            testCase.verifyEqual(ismember(24,samples),true);
            testCase.verifyEqual(ismember(48,samples),true);
        end

        function MIMOfixedSamplePoints(testCase,sysMIMOn40)
        % Use fixed sample points
            opts.maxIter = 10;
            opts.adaptiveSampling = false;
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.method, @sobmor);
            testCase.verifyEqual(length(redSys.parameters.samplePoints), 807);
        end

        function MIMOgammaSequence(testCase,sysMIMOn40)
        % Use fixed sample points
            opts.maxIter = 10;
            opts.gammaSequence = [5,4,3,2,1];
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);
            testCase.verifyEqual(redSys.method, @sobmor);
            testCase.verifyEqual(length(redSys.parameters.gammaSequence), 5);
        end

        function MIMOGransoOptions(testCase,sysMIMOn40)
        % Use fixed sample points
            opts.maxIter = 10;
            opts.granso.step_tol = 1e-10;
            opts.granso.maxit = 50;
            opts.granso.print_level = 3;
            redSys = sobmor(sysMIMOn40,6,opts);
            testCase.verifyEqual(redSys.dim, 6);

            % sobmor overwrites max iterations and print level input to granso
            testCase.verifyEqual(redSys.parameters.granso.maxit, 10);
            testCase.verifyEqual(redSys.parameters.granso.step_tol, 1e-10);
            testCase.verifyEqual(redSys.parameters.granso.print_level, 0);
        end

    end

    methods (Test, TestTags = {'Functionality'})
        % Simulate more realistic scenarios with larger original models
        % Compare results with other functions (irkaPH)
        function accuracy(testCase, sysSISOn300)
        % make sure reduced system approximates the original
            redSys_SISO = sobmor(sysSISOn300,10);
            testCase.verifyLessThanOrEqual(norm(redSys_SISO-sysSISOn300,Inf)/norm(sysSISOn300,Inf), 1e-2);
            
        % Compare results of sobmor to those of irkaPH
            opts.summary = false;
            redSys_irka = irkaPH(sysSISOn300,10,opts);
            err_irka = norm(redSys_irka-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            err_sobmor = norm(redSys_SISO-sysSISOn300,Inf)/norm(sysSISOn300,Inf);
            testCase.verifyLessThan(err_sobmor, err_irka);
        end

        function sobmorVsIrkaPHMIMO(testCase, sysMIMOn300)
            % Compare results of sobmor to those of irkaPH (MIMO case)
            opts.constAlg = 'phs';
            opts.tol = 1e-10;
            opts.summary = false;
            redSys_irka = irkaPH(sysMIMOn300,10,opts);
            redSys_MIMO = sobmor(setup_MassSpringDamperSystem(300,2,1,1,'MIMO'),10);
            err_irka = norm(redSys_irka-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            err_sobmor = norm(redSys_MIMO-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            testCase.verifyLessThanOrEqual(err_sobmor, err_irka);
            
            % Compare results of sobmor to those of Effort Constraint Method (MIMO case)            
            redSys_ecm = effortConstraint(sysMIMOn300,10);
            err_ecm = norm(redSys_ecm-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            err_sobmor = norm(redSys_MIMO-sysMIMOn300,Inf)/norm(sysMIMOn300,Inf);
            testCase.verifyLessThanOrEqual(err_sobmor, err_ecm);
        end
    end
end