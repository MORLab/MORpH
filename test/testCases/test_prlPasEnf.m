classdef test_prlPasEnf < matlab.unittest.TestCase
    % TEST test_prlPasEnf - Tests the functionalities of dominant spectral zero
    % method
    %
    % Press <F5> or enter "runtests("test_prlPasEnf")" to run this testscript
    %
    % Description:
    %   This script tests the function 'prlPasEnf' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Functionality: Check passivity of perturbed models for randomized
    %                        non-passive systems
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
        OptsSP = struct('plot',false,'outputVolume','worst');
        sys = {sss(setup_MassSpringDamperSystem(100,4,4,1,'MIMO'))};
        sysr4 = {irka(sss(setup_MassSpringDamperSystem(100,4,4,1,'MIMO')),4,struct('maxiter',50,'tol',1e-3,'stopCrit','combAny'))};
    end

    properties

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
    %% Test in and outputs
    methods (Test, TestTags = {'Input','Output'})
        function standardInput(test,sysr4)
            % Very basic input: system
            sysp = prlPasEnf(sysr4);
            test.verifyEqual(size(sysp.A,1), 4);
            test.verifyEqual(sysp.redParam.method, 'prlPasEnf');
            test.verifyTrue(isa(sysp,'ssRed'));
            test.verifyTrue(samPassive(sysp,test.OptsSP));
        end

        % Input
        function standardInputWithOpts(testCase,sysr4)
            Opts.makePH = true;
            Opts.tol = 1e-6;

            sysp = prlPasEnf(sysr4,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),Opts.tol);
            testCase.verifyEqual(size(sysp.Q,1), 4);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function cvx(testCase,sysr4)
            Opts.solver = 'cvx';
            Opts.makePH = true;
            Opts.tol = 1e-6;

            sysp = prlPasEnf(sysr4,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),Opts.tol);
            testCase.verifyEqual(size(sysp.Q,1), 4);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end
        
        function yalmip(testCase,sysr4)            
            Opts.solver = 'yalmip';
            Opts.makePH = true;
            
            sysp = prlPasEnf(sysr4,Opts);
            
            testCase.verifyEqual(size(sysp.Q,1), 4);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end
    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function Passivity4(testCase,sysr4)
            Opts.makePH = true;

            sysp = prlPasEnf(sysr4,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),1e-8);
            testCase.verifyEqual(size(sysp.Q,1), 4);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function Passivity10(testCase,sys)
            sysr = irka(sys,10);

            Opts.makePH = true;

            sysp = prlPasEnf(sysr,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),1e-8);
            testCase.verifyEqual(size(sysp.Q,1), 10);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function Passivity16(testCase,sys)
            sysr = irka(sys,16);

            Opts.makePH = true;

            sysp = prlPasEnf(sysr,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),1e-8);
            testCase.verifyEqual(size(sysp.Q,1), 16);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function Passivity8(testCase)
            options.ft = true;
            options.P = 'rand';
            rng(12);
            sysnp = setup_RandomSystem(8,1,-1e-1,options);

            Opts.makePH = true;

            sysp = prlPasEnf(sysnp,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),1e-8);
            testCase.verifyEqual(size(sysp.Q,1), 8);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function Passivity12(testCase)
            options.ft = true;
            options.P = 'rand';
            rng(8);
            sysnp = setup_RandomSystem(12,1,-1e-1,options);

            Opts.makePH = true;

            sysp = prlPasEnf(sysnp,Opts);

            d = eig(sysp.parameters.X);
            testCase.verifyGreaterThan(min(d),1e-8);
            testCase.verifyEqual(size(sysp.Q,1), 12);
            testCase.verifyTrue(isa(sysp,'phsRed'));
        end

        function Passivity20(testCase)
            options.ft = true;
            options.P = 'rand';
            rng(1);
            sysr = setup_RandomSystem(20,1,-1e-3,options);

            sysp = prlPasEnf(sysr);

            testCase.verifyEqual(size(sysp.A,1), 20);
            testCase.verifyTrue(samPassive(sysp,struct('plot',false)));
        end
    end
end