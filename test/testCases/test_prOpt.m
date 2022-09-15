classdef test_prOpt < matlab.unittest.TestCase
    % TEST test_prOpt - Tests the functionalities of <module name>
    %
    % Press <F5> or enter "runtests("test_prOpt")" to run this testscript
    %
    % Description:
    %   This script tests the function 'prOpt' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Output: Check if each additional output can be obtained
    %       - Options: Check if additional options work correctly
    %
    %-----------------------------------------------------------------------
    % This file is part of
    %
    % <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
    % interconnect and reduce large-scale port-Hamiltonian models
    %
    % Authors:      Tim Moser
    % E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
    % Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
    % Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
    %-----------------------------------------------------------------------

    %% Properties

    properties (TestParameter) % test parameters
        sys1 = {setup_LadderNetworkSystem(20,0.1,0.1,1)};
        sysM = {setup_MassSpringDamperSystem(10,4,4,1,'MIMO')};
        sysDAE1 = {setup_RandomStaircasePHDAE([0,10,4,0],2)};
        sysDAE12 = {setup_RandomStaircasePHDAE([1,10,4,1],1)};
        sysDAE2 = {setup_RandomStaircasePHDAE([1,10,0,1],1)};
    end

    properties % other properties (not used as test inputs)
        tol = 1e-10;
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
    % You may also combine test tags for some functions

    methods (Test, TestTags = {'Input'})
        % Check input combinations

        % Test standardInput (original system and reduced-order)
        function standardInputs(testCase,sys1)
            Opts.maxIter = 1;
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.dim,4);
        end

        % Test handover of initial model
        function initialModel(testCase,sysDAE12)
            redSys0 = pHMOR_DAE_Wrapper(@irkaPH,sysDAE12,2);
            Opts.maxIter = 0;
            redSys = prOpt(sysDAE12,2,redSys0,Opts);
            testCase.verifyLessThan(norm(redSys0-redSys), testCase.tol);

            Opts.optParams = struct('Er',false,'Jr',false,'Rr',false,'Qr',false,'Gr',false,'Pr',false);
            redSys = prOpt(sysDAE12,2,redSys0,Opts);
            testCase.verifyLessThan(norm(redSys0-redSys), testCase.tol);
        end

    end

    % Check if (optional) output is correct
    methods (Test, TestTags = {'Output'})

        % Test if reduced model has port-Hamiltonian form (if input
        % Validation is not performed for some reason)
        function isPH(testCase,sys1)
            Opts.maxIter = 10;
            redSys = prOpt(sys1,6,Opts);
            isPH = redSys.inputValidation(redSys);
            testCase.verifyEqual(isPH,true);
        end

    end

    % Check if options are applied correctly
    methods (Test, TestTags = {'Options'})

        % Check different optimization parameter configurations
        function optParams(testCase,sys1)
            Opts.maxIter = 3;
            % default
            Opts.optParams = struct('Er',false,'Jr',true,'Rr',true,'Qr',true,'Gr',true,'Pr',true);
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyLessThan(norm(redSys.E-eye(4)), testCase.tol);
            % E,J,R,G
            Opts.optParams = struct('Er',false,'Jr',true,'Rr',true,'Qr',false,'Gr',true,'Pr',false);
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyLessThan(norm(redSys.Q-eye(4)), testCase.tol);
            testCase.verifyLessThan(norm(redSys.P-zeros(4,1)), testCase.tol);
            % J,R
            Opts.optParams = struct('Er',false,'Jr',true,'Rr',true,'Qr',false,'Gr',false,'Pr',false);
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyLessThan(norm(redSys.E-eye(4)), testCase.tol);
            testCase.verifyLessThan(norm(redSys.Q-eye(4)), testCase.tol);
            testCase.verifyLessThan(norm(redSys.G-zeros(4,1)), testCase.tol);
            testCase.verifyLessThan(norm(redSys.P-zeros(4,1)), testCase.tol);
        end

        % Check different algorithms for EVP differentiation
        function eigDerivative(testCase,sys1)
            % Test options for Lyapunov solver
            Opts.maxIter = 10;
            Opts.eigDerivative = 'magnus';
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.parameters.eigDerivative,'magnus');
            Opts.eigDerivative = 'rudisill';
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.parameters.eigDerivative,'rudisill');
            Opts.eigDerivative = 'murthy';
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.parameters.eigDerivative,'murthy');
        end

        function maxIter(testCase,sys1)
            Opts.maxIter = 3;
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyEqual(3,redSys.info.convergence(end).iter);
        end

        function tolGradNorm(testCase,sys1)
            Opts.tol = 1e-2;
            redSys = prOpt(sys1,4,Opts);
            testCase.verifyLessThanOrEqual(redSys.info.convergence(end).gradnorm, 1e-2);
            testCase.verifyGreaterThan(redSys.info.convergence(end).gradnorm, 1e-6);
        end

        % Check different solvers
        function solver(testCase,sysM)
            Opts.maxIter = 3;
            Opts.solver = 'Manopt';
            redSys = prOpt(sysM,4,Opts);
            testCase.verifyEqual(redSys.parameters.solver,'Manopt');
            Opts.solver = 'GRANSO';
            redSys = prOpt(sysM,4,Opts);
            testCase.verifyEqual(redSys.parameters.solver,'GRANSO');
            Opts.solver = 'MATLAB';
            redSys = prOpt(sysM,4,Opts);
            testCase.verifyEqual(redSys.parameters.solver,'MATLAB');
        end

    end

    % Check if functionality is given (e.g. accuracy of results)
    methods (Test, TestTags = {'Functionality'})

        % Test the Riemannian gradient computation
        function checkGradient(testCase,sys1)
            Opts.maxIter = 1;
            Opts.manopt.checkGradient = true;
            redSys = prOpt(sys1,8,Opts);
            userCheck = input('Do you accept the slope of the gradient? y/n [y]: ','s');
            if isempty(userCheck)
                userCheck = 'n';
            end
            testCase.verifyEqual(userCheck,'y');
            close;
        end

        % Test the Riemannian Hessian computation
        function checkHessian(testCase,sys1)
            Opts.maxIter = 1;
            Opts.manopt.checkHessian = true;
            redSys = prOpt(sys1,4,Opts);
            userCheck = input('Do you accept the slope of the Hessian? y/n [y]: ','s');
            if isempty(userCheck)
                userCheck = 'n';
            end
            testCase.verifyEqual(userCheck,'y');
            close;
        end

    end

    methods (Test, TestTags = {'Benchmarks'})
        % Use benchmarks if available

    end

    %% Supporting functions (e.g. for reuse)

    methods

    end

end