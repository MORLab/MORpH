classdef test_ss2phs < matlab.unittest.TestCase
    % TEST test_ss2phs - Tests the functionalities of positive real and mixed
    % gramian balancing
    %
    % Press <F5> or enter "runtests("test_ss2phs")" to run this testscript
    %
    % Description:
    %   This script tests the function 'ss2phs' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Functionality: Check if reduced models are transformed to pH
    %                        systems
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
        sys100 = {sss(setup_MassSpringDamperSystem(100,4,4,1))};      % MSD n=40
        sysr4 = {sss(irkaPH(setup_MassSpringDamperSystem(100,4,4,1),4))};
        sys100MIMO = {sss(setup_MassSpringDamperSystem(100,4,4,1,'MIMO'))};      % MSD n=100
        sysr4MIMO = {sss(irkaPH(setup_MassSpringDamperSystem(100,4,4,1,'MIMO'),4))};
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
        function standardInput(test,sysr4)
            % Very basic input: system only
            sysr = ss2phs(sysr4);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.method, @ss2phs);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function standardInputMIMO(test,sysr4MIMO)
            % Very basic input: system only
            sysr = ss2phs(sysr4MIMO);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.method, @ss2phs);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function standardWithOpts(test,sysr4)
            Opts.method = @irkaPH;
            Opts.prl = 'icare';
            Opts.representation = 'standard';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.method, @irkaPH);
            test.verifyEqual(sysr.parameters.prl, 'icare');
            test.verifyEqual(sysr.parameters.representation, 'standard');
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        % Representations
        function implicit(test,sysr4)
            Opts.representation = 'implicit';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.parameters.representation, 'implicit');
            test.verifyTrue(isa(sysr,'phsRed'));
            test.verifyEqual(sysr.Q, eye(size(sysr.Q)));
            test.verifyNotEqual(sysr.E, eye(size(sysr.E)));
        end

        function standard(test,sysr4)
            Opts.representation = 'standard';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.parameters.representation, 'standard');
            test.verifyTrue(isa(sysr,'phsRed'));
            test.verifyNotEqual(sysr.Q, eye(size(sysr.Q)));
            test.verifyEqual(sysr.E, eye(size(sysr.E)));
        end

        function scaled(test,sysr4)
            Opts.representation = 'scaled';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyEqual(sysr.parameters.representation, 'scaled');
            test.verifyTrue(isa(sysr,'phsRed'));
            test.verifyEqual(sysr.Q, eye(size(sysr.Q)));
            test.verifyEqual(sysr.E, eye(size(sysr.E)));
        end

        % Solver
        function cvx(test,sysr4)
            Opts.prl = 'cvx';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function icare(test,sysr4)
            Opts.prl = 'icare';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function mess(test,sys100)
            Opts.phs.inputTolerance = 1e-5;
            Opts.prl = 'mmess';
            sysr = ss2phs(sys100,Opts);
            test.verifyEqual(size(sysr.J,1), 100);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function yalmip(test,sysr4)
            Opts.prl = 'yalmip';
            sysr = ss2phs(sysr4,Opts);
            test.verifyEqual(size(sysr.J,1), 4);
            test.verifyTrue(isa(sysr,'phsRed'));
        end
    end

    %% Check if reduced system is pH
    methods (Test, TestTags = {'Functionality'})
        % SISO
        function r6(test,sys100)
            sysr = irka(sys100,6,struct('maxiter',200,'tol',1e-6,'stopCrit','s0'));
            sysr = ss2phs(sysr);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function r20(test,sys100)
            sysr = irka(sys100,20,struct('maxiter',200,'tol',1e-6,'stopCrit','s0'));
            sysr = ss2phs(sysr);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        % MIMO
        function r6MIMO(test,sys100)
            sysr = prbt(sys100,6);
            sysr = ss2phs(sysr);
            test.verifyTrue(isa(sysr,'phsRed'));
        end

        function r20MIMO(test,sys100MIMO)
            sysr = irka(sys100MIMO,20,struct('maxiter',200,'tol',1e-6,'stopCrit','s0'));
            sysr = ss2phs(sysr);
            test.verifyTrue(isa(sysr,'phsRed'));
        end
    end
end