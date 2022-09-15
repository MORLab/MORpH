classdef test_hamPasEnf < matlab.unittest.TestCase
    % TEST test_hamPasEnf - Tests the functionalities of dominant spectral zero
    % method
    %
    % Press <F5> or enter "runtests("test_hamPasEnf")" to run this testscript
    %
    % Description:
    %   This script tests the function 'hamPasEnf' of the MORpH-toolbox. The
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
        function standardInput(test)
            % Very basic input: system
            rng('default');
            sys = setup_RandomSystem(4,1,-1e-3);

            sysp = hamPasEnf(sys);

            test.verifyEqual(size(sysp.A,1), 4);
            test.verifyEqual(sysp.redParam.method, 'hamPasEnf');
            test.verifyTrue(isa(sysp,'ssRed'));
            test.verifyTrue(samPassive(sysp,test.OptsSP));
        end

        % Input
        function standardInputWithOpts(testCase)
            rng('default');
            sys = setup_RandomSystem(4,1,-1e-3);

            Opts.makePH = true;
            Opts.maxIter = 500;
            Opts.verbose = 1;
            Opts.alpha = 0.4;
            Opts.ss2phs.phs.inputTolerance = 1e-7;
            Opts.ss2phs.prl = 'icare';
            sysp = hamPasEnf(sys,Opts);

            testCase.verifyEqual(size(sysp.Q,1), 4);
            testCase.verifyTrue(isa(sysp,'phsRed'));
            testCase.verifyEqual(sysp.parameters.maxIter, 500);
            testCase.verifyEqual(sysp.parameters.alpha, 0.4);
        end

    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function Passivity4(testCase)
            options.ft = true;
            options.P = 'rand';
            rng('default');
            sys = setup_RandomSystem(4,1,-1e-3,options);

            sysp = hamPasEnf(sys);

            testCase.verifyEqual(size(sysp.A,1), 4);
            testCase.verifyTrue(samPassive(sysp,testCase.OptsSP));
        end

        function Passivity8(testCase)
            options.ft = true;
            options.P = 'rand';
            rng('default');
            sys = setup_RandomSystem(8,1,-1e-3,options);

            sysp = hamPasEnf(sys);

            testCase.verifyEqual(size(sysp.A,1), 8);
            testCase.verifyTrue(samPassive(sysp,testCase.OptsSP));
        end

        function Passivity12(testCase)
            options.ft = true;
            options.P = 'rand';
            rng('default');
            sys = setup_RandomSystem(12,1,-1e-3,options);

            sysp = hamPasEnf(sys);

            testCase.verifyEqual(size(sysp.A,1), 12);
            testCase.verifyTrue(samPassive(sysp,testCase.OptsSP));
        end

        function Passivity16(testCase)
            options.ft = true;
            options.P = 'rand';
            rng('default');
            sys = setup_RandomSystem(16,1,-1e-3,options);

            sysp = hamPasEnf(sys);

            testCase.verifyEqual(size(sysp.A,1), 16);
            testCase.verifyTrue(samPassive(sysp,testCase.OptsSP));
        end


        function Passivity20(testCase)
            options.ft = true;
            options.P = 'rand';
            rng('default');
            sys = setup_RandomSystem(20,1,-1e-3,options);

            sysp = hamPasEnf(sys);

            testCase.verifyEqual(size(sysp.A,1), 20);
            testCase.verifyTrue(samPassive(sysp,testCase.OptsSP));
        end
    end
end