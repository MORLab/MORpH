classdef test_lyapOpt < matlab.unittest.TestCase
    % TEST test_lyapOpt - Tests the functionalities of <module name>
    %
    % Press <F5> or enter "runtests("test_lyapOpt")" to run this testscript
    %
    % Description:
    %   This script tests the function 'lyapOpt' of the MORpH-toolbox. The
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
        sys2 = {setup_MassSpringDamperSystem(10,4,4,1,'MIMO')};
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

    % Check input combinations
    methods (Test, TestTags = {'Input'})

        % Test standardInput (original system and reduced-order)
        function standardInputs(testCase,sys1)
            Opts.manopt.maxiter = 10;
            redSys = lyapOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.dim,4);
        end

        % Test handover of initial model
        function initialModel(testCase,sys2)
            redSys0 = arnoldiPH(sys2,zeros(1,4));
            Opts.manopt.maxiter = 0;
            redSys = lyapOpt(sys2,8,redSys0,Opts);
            testCase.verifyLessThan(norm(redSys0-redSys), testCase.tol);
        end
    end

    methods (Test, TestTags = {'Output'})

        % Test if reduced model has port-Hamiltonian form (if input
        % Validation is not performed for some reason)
        function isPH(testCase,sys2)
            Opts.manopt.maxiter = 20;
            redSys = lyapOpt(sys2,2,Opts);
            isPH = redSys.inputValidation(redSys);
            testCase.verifyEqual(isPH,true);
        end

    end

    methods (Test, TestTags = {'Options'})
        % Check if options are applied correctly

        function inp_lyapSolver(testCase,sys1)
            % Test options for Lyapunov solver
            Opts.manopt.maxiter = 10;
            Opts.lyapSolver = 'mess';
            redSys = lyapOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.parameters.lyapSolver,'mess');
            Opts.lyapSolver = 'lyap';
            redSys = lyapOpt(sys1,4,Opts);
            testCase.verifyEqual(redSys.parameters.lyapSolver,'lyap');
        end

        function maxIter(testCase,sys1)
            Opts.manopt.maxiter = 3;
            redSys = lyapOpt(sys1,4,Opts);
            testCase.verifyEqual(3,redSys.info(end).iter);

        end
        function tolGradNorm(testCase,sys2)
            Opts.manopt.tolgradnorm = 1e-2;
            redSys = lyapOpt(sys2,4,Opts);
            testCase.verifyLessThanOrEqual(redSys.info(end).gradnorm, 1e-2);
            testCase.verifyGreaterThan(redSys.info(end).gradnorm, 1e-6);

        end
    end

    % Check if functionality is given (e.g. accuracy of results)
    methods (Test, TestTags = {'Functionality'})

        % Test the Riemannian gradient computation
        function func_checkGradient_lyap(testCase,sys1)
            Opts.manopt.maxiter = 1;
            Opts.checkGradient = true;
            Opts.lyapSolver = 'lyap';
            lyapOpt(sys1,8,Opts);
            userCheck = input('Do you accept the slope of the gradient? y/n [y]: ','s');
            if isempty(userCheck)
                userCheck = 'n';
            end
            testCase.verifyEqual(userCheck,'y');
            close;
        end

        % Test the Riemannian Hessian computation
        function func_checkHessian_lyap(testCase,sys2)
            Opts.manopt.maxiter = 1;
            Opts.checkHessian = true;
            Opts.lyapSolver = 'lyap';
            lyapOpt(sys2,4,Opts);
            userCheck = input('Do you accept the slope of the Hessian? y/n [y]: ','s');
            if isempty(userCheck)
                userCheck = 'n';
            end
            testCase.verifyEqual(userCheck,'y');
            close;
        end

        % Test the Riemannian gradient computation
        function func_checkGradient_mess(testCase,sys1)
            Opts.manopt.maxiter = 1;
            Opts.checkGradient = true;
            Opts.lyapSolver = 'mess';
            lyapOpt(sys1,8,Opts);
            userCheck = input('Do you accept the slope of the gradient? y/n [y]: ','s');
            if isempty(userCheck)
                userCheck = 'n';
            end
            testCase.verifyEqual(userCheck,'y');
            close;
        end

        % Test the Riemannian Hessian computation
        function func_checkHessian_mess(testCase,sys2)
            Opts.manopt.maxiter = 1;
            Opts.checkHessian = true;
            Opts.lyapSolver = 'mess';
            lyapOpt(sys2,4,Opts);
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