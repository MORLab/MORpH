classdef test_dszm < matlab.unittest.TestCase
    % TEST test_dszm - Tests the functionalities of dominant spectral zero
    % method
    %
    % Press <F5> or enter "runtests("test_dszm")" to run this testscript
    %
    % Description:
    %   This script tests the function 'dszm' of the MORpH-toolbox. The
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
        sys100noFT = {setup_RandomPassiveSystem(100,1,struct('P','zero','ft',false))};
        sys100FT = {setup_RandomPassiveSystem(100,1,struct('P','rand','ft',true))};
    end

    properties

    end

    %% Setup / teardown of test environment

    methods (TestClassSetup) % once before all tests
        function set_rng(testCase)
            rng(1)
        end
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
        function standardInput(test,sys40)
            % Very basic input: system and reduced order only
            sysr = dszm(sys40,4);
            test.verifyLessThanOrEqual(size(sysr.A,1), 5);
            test.verifyEqual(sysr.redParam.method, 'dszm');
            test.verifyTrue(isa(sysr,'ssRed'));
        end

        % Input
        function standardInputWithOpts(testCase,sys40)
            Opts.tol = 1e-10;
            Opts.explicit = true;
            Opts.makePH = true;
            Opts.checkPassivity = false;
            Opts.projectiveMOR.trans = 'H';
            Opts.arnoldi.orth = 'dgks';
            Opts.samPassive.tol = 1e-6;
            Opts.ss2phs.representation = 'standard';

            sysr = dszm(sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.parameters.tol, 1e-10);
            testCase.verifyTrue(sysr.parameters.explicit);
            testCase.verifyFalse(sysr.parameters.checkPassivity);
            testCase.verifyMatches(sysr.parameters.projectiveMOR.trans,'H');
            testCase.verifyMatches(sysr.parameters.arnoldi.orth, 'dgks');
            testCase.verifyEqual(sysr.parameters.samPassive.tol, 1e-6);
            testCase.verifyEqual(sysr.E, eye(size(sysr.E)));
        end

    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function Passivity4(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,4,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity6(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,6,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity8(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,8,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity10(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,10,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity12(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,12,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity14(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,14,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity16(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,16,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity18(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,18,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity20(testCase,sys100FT)
            Options.makePH = true;
            sysr = dszm(sys100FT,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        % No feed-through
        function Passivity4_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,4,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity6_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,6,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity8_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,8,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity10_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,10,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity12_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,12,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity14_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,14,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity16_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,16,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity18_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,18,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end

        function Passivity20_no_FT(testCase,sys100noFT)
            Options.makePH = true;
            sysr = dszm(sys100noFT,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
        end
    end
end