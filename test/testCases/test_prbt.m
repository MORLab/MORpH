classdef test_prbt < matlab.unittest.TestCase
    % TEST test_prbt - Tests the functionalities of positive real and mixed
    % Gramian balancing
    %
    % Press <F5> or enter "runtests("test_prbt")" to run this testscript
    %
    % Description:
    %   This script tests the function 'prbt' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Functionality: Check passivity of reduced models for randomized
    %                        passive systems
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
        sys40 = {sss(setup_MassSpringDamperSystem(40,4,4,1))};      % MSD n=40
        sys100FT = {setup_RandomPassiveSystem(100,1,struct('P','rand','ft',true))};
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

            if ~(exist('setup_RandomPassiveSystem.m','file') == 2)
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
        function standardInput(test,sys40)
            % Very basic input: system and reduced order only
            sysr = prbt(sys40,4);
            test.verifyEqual(size(sysr.A,1), 4);
            test.verifyEqual(sysr.redParam.method, 'prbt');
            test.verifyTrue(isa(sysr,'ssRed'));
        end

        % Input
        function standardInputWithOpts(testCase,sys40)
            Opts.useCholFactors = false;
            Opts.makePH = true;
            Opts.ss2phs.representation = 'standard';
            Opts.samPassive.tol = 1e-6;
            Opts.tol = 1e-10;
            Opts.checkPassivity = false;
            Opts.mixedGramian = 'standard';
            Opts.lyap = 'lyapchol';
            Opts.are = 'icare';

            sysr = prbt(sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.parameters.tol, 1e-10);
            testCase.verifyFalse(sysr.parameters.checkPassivity);
            testCase.verifyFalse(sysr.parameters.useCholFactors);
            testCase.verifyMatches(sysr.parameters.mixedGramian,'standard');
            testCase.verifyMatches(sysr.parameters.lyap,'lyapchol');
            testCase.verifyEqual(sysr.parameters.samPassive.tol, 1e-6);
            testCase.verifyEqual(sysr.E, eye(size(sysr.E)));
        end

        function mixedControlMMESS(testCase,sys40)
            Opts.mixedGramian = 'mixedControl';
            Opts.are = 'mmess';
            Opts.lyap = 'mmess';
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function mixedControl(testCase,sys40)
            Opts.mixedGramian = 'mixedControl';
            Opts.are = 'icare';
            Opts.lyap = 'lyapchol';
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function mixedObserveMMESS(testCase,sys40)
            Opts.mixedGramian = 'mixedObserve';
            Opts.are = 'mmess';
            Opts.lyap = 'mmess';
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function mixedObserve(testCase,sys40)
            Opts.mixedGramian = 'mixedObserve';
            Opts.are = 'icare';
            Opts.lyap = 'lyapchol';
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function User(testCase,sys40)
            Opts.truncation = 'userIn';
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function TruncTol(testCase,sys40)
            Opts.truncation = 'truncTol';
            Opts.truncTol = 1e-3;
            Opts.makePH = true;
            sysr = prbt(sys40,4,Opts);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end
    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function Passivity4(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,4,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity6(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,6,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity8(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,8,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity10(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,10,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity12(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,12,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity14(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,14,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity16(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,16,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity18(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,18,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity20(testCase,sys100FT)
            Options.makePH = true;
            sysr = prbt(sys100FT,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity20MixedControl(testCase,sys100FT)
            Options.mixedGramian = 'mixedControl';
            Options.makePH = true;
            sysr = prbt(sys100FT,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity20MixedObserve(testCase,sys100FT)
            Options.mixedGramian = 'mixedObserve';
            Options.makePH = true;
            sysr = prbt(sys100FT,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end
    end
end