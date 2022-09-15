classdef test_balPH < matlab.unittest.TestCase
    % TEST test_balPH - Tests the functionalities of dominant spectral zero
    % method
    %
    % Press <F5> or enter "runtests("test_balPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'balPH' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Functionality: Check if reduced models are port-Hamiltonian
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
        sys40 = {setup_MassSpringDamperSystem(40,4,4,1)};      % MSD n=40
        sys100_MSD = {setup_MassSpringDamperSystem(100,4,4,1)};
        sys100_LN = {setup_LadderNetworkSystem(100,.1,.1,3)};
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
        function standardInput(test,sys40)
            % Very basic input: system and reduced order only
            sysr = balPH(sys40);
            test.verifyEqual(sysr.method, @balPH);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function reduced_order(test,sys40)
            % Very basic input: system and reduced order only
            sysr = balPH(sys40,4);
            test.verifyEqual(size(sysr.R,1), 4);
            test.verifyEqual(sysr.method, @balPH);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        % Input
        function standardInputWithOpts(testCase,sys40)
            Opts.truncTol = 0.1;
            Opts.lyap = 'lyapchol';

            sysr = balPH(sys40,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.R,1), 14);
            testCase.verifyEqual(sysr.parameters.truncTol, 0.1);
            testCase.verifyEqual(sysr.parameters.lyap, 'lyapchol');
        end

    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function PHS_MSD_4(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,4,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),4);
        end

        function PHS_MSD_6(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,6,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),6);
        end

        function PHS_MSD_8(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,8,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),8);
        end

        function PHS_MSD_10(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,10,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),10);
        end

        function PHS_MSD_12(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,12,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),12);
        end

        function PHS_MSD_14(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,14,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),14);
        end

        function PHS_MSD_16(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,16,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),16);
        end

        function PHS_MSD_18(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,18,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),18);
        end

        function PHS_MSD_20(testCase,sys100_MSD)
            Options.makePH = true;
            sysr = balPH(sys100_MSD,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),20);
        end

        function PHS_LN_4(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,4,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),4);
        end

        function PHS_LN_6(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,6,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),6);
        end

        function PHS_LN_8(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,8,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),8);
        end

        function PHS_LN_10(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,10,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),10);
        end

        function PHS_LN_12(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,12,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),12);
        end

        function PHS_LN_14(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,14,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),14);
        end

        function PHS_LN_16(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,16,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),16);
        end

        function PHS_LN_18(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,18,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),18);
        end

        function PHS_LN_20(testCase,sys100_LN)
            Options.makePH = true;
            sysr = balPH(sys100_LN,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),20);
        end

        function PHS_FT_m1(testCase)
            Options.makePH = true;
            rng(1)
            [~,sys]=setup_RandomPassiveSystem(100,1,struct('P','rand','ft',true));
            sysr = balPH(sys,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),20);
        end

        function PHS_FT_m2(testCase)
            Options.makePH = true;
            rng(1)
            [~,sys]=setup_RandomPassiveSystem(100,2,struct('P','rand','ft',true));
            sysr = balPH(sys,20,Options);
            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(size(sysr.J,1),20);
        end
    end
end