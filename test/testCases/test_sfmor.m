classdef test_sfmor < matlab.unittest.TestCase
    % TEST test_sfmor - Tests the functionalities of positive real and mixed
    % gramian balancing
    %
    % Press <F5> or enter "runtests("test_sfmor")" to run this testscript
    %
    % Description:
    %   This script tests the function 'sfmor' of the MORpH-toolbox. The
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
            % Very basic input: system only
            [sysr,~,rfac] = sfmor(@irka,sys40,10);
            test.verifyEqual(size(sysr.A,1), 10);
            test.verifyEqual(sysr.redParam.method, 'sfmor');
            test.verifyTrue(isa(sysr,'ssRed'));
            test.verifyEqual(rfac.redParam.method,'irka');
        end

        % Input
        function standardInputWithOpts(testCase,sys40)
            Opts.useCholFactors = false;
            Opts.makePH = true;
            Opts.ss2phs.representation = 'standard';
            Opts.samPassive.tol = 1e-6;
            Opts.tol = 1e-10;
            Opts.checkPassivity = false;
            Opts.lyap = 'lyap';
            Opts.are = 'icare';

            [sysr,~,rfac] = sfmor(@tbr,sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.parameters.tol, 1e-10);
            testCase.verifyFalse(sysr.parameters.checkPassivity);
            testCase.verifyFalse(sysr.parameters.useCholFactors);
            testCase.verifyMatches(sysr.parameters.lyap,'lyap');
            testCase.verifyEqual(sysr.parameters.samPassive.tol, 1e-6);
            testCase.verifyEqual(sysr.E, eye(size(sysr.E)));
            testCase.verifyEqual(rfac.redParam.method,'tbr');
        end

        function MMESS(testCase,sys40)
            Opts.lyap = 'mmess';
            Opts.are = 'mmess';
            Opts.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyMatches(sysr.parameters.lyap,'mmess');
            testCase.verifyMatches(sysr.parameters.are,'mmess');
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function LyapIcare(testCase,sys40)
            Opts.lyap = 'lyap';
            Opts.are = 'icare';

            [sysr,~,rfac] = sfmor(@irka,sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'ssRed'));
            testCase.verifyMatches(sysr.redParam.params.lyap,'lyap');
            testCase.verifyMatches(sysr.redParam.params.are,'icare');
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function AdditionalInputIRKA(testCase,sys40)
            Opts.irka.tol=1e-4;
            Opts.irka.stopCrit = 's0';
            Opts.irka.maxiter = 200;
            Opts.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys40,4,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(rfac.redParam.method,'irka');
            testCase.verifyEqual(rfac.redParam.params.tol,1e-4);
            testCase.verifyEqual(rfac.redParam.params.stopCrit,'s0');
            testCase.verifyEqual(rfac.redParam.params.maxiter,200);
        end

        function UserInTbr(testCase,sys40)
            Opts.makePH = true;

            [sysr,~,rfac] = sfmor(@tbr,sys40,Opts);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(rfac.redParam.method,'tbr');
        end
    end

    %% Check if reduced system is passive and pH
    methods (Test, TestTags = {'Functionality'})
        function Passivity4(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,4,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,4);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity6(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,6,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,6);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity8(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,8,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,8);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity10(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,10,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,10);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity12(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,12,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,12);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity14(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,14,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,14);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity16(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,16,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,16);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity18(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,18,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,18);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity20(testCase,sys100FT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100FT,20,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,20);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        % No Feed-Through
        function Passivity4_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,4,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,4);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity6_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,6,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,6);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity8_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,8,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,8);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity10_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,10,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,10);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity12_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,12,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,12);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity14_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,14,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,14);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity16_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,16,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,16);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity18_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,18,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,18);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end

        function Passivity20_noFT(testCase,sys100noFT)
            Options.makePH = true;

            [sysr,~,rfac] = sfmor(@irka,sys100noFT,20,Options);

            testCase.verifyTrue(isa(sysr,'phsRed'));
            testCase.verifyEqual(sysr.dim,20);
            testCase.verifyEqual(rfac.redParam.method,'irka');
        end
    end
end