classdef test_trksmPH < matlab.unittest.TestCase
    % TEST test_trksmPH - Tests the functionalities of trksmPH
    %
    % Press <F5> or enter "runtests("test_trksmPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'trksmPH' of the MORpH-toolbox. The
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
    % Authors:      Tim Moser
    % E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
    % Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
    % Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
    %-----------------------------------------------------------------------

    %% Properties
    properties (TestParameter)
        sys = {setup_MassSpringDamperSystem(40,2,1,1)};   % standard MSD system for most of the tests
        sysMIMO = {setup_MassSpringDamperSystem(40,2,1,1,'MIMO')}; % MIMO MSD system for MIMO testing
        shifts = {logspace(-5,5,5)}; % standard shifts
    end

    properties
        tolerance = 1e-10;
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


    %% Tests
    methods (Test, TestTags = {'Input'})
        % check input combinations

        function standardInput(testCase,sys)
            % Test standard input (system + reduced order)

            redSys = trksmPH(sys,5);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function standardInputWithOpts(testCase,sys)
            % Test standard input with simple options

            Opts.tol = 1e-2;
            Opts.minSamplePoints = 100;
            Opts.isReal = true;
            redSys = trksmPH(sys,5,Opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.tol, Opts.tol);
            testCase.verifyEqual(redSys.parameters.minSamplePoints, Opts.minSamplePoints);
            testCase.verifyEqual(redSys.parameters.isReal, Opts.isReal);
        end

        function initialShifts(testCase,sys)
            % Supply initial shifts

            redSys = trksmPH(sys, 5, 1e-5,1e-3);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyTrue(any(redSys.info.shifts==1e-5));
            testCase.verifyEqual(redSys.parameters.sMin, 1e-5);
            testCase.verifyEqual(redSys.parameters.sMax, 1e-3);
        end

        function MIMO(testCase,sysMIMO)
            % Test MIMO system reduction
            redSys = trksmPH(sysMIMO,10);
            testCase.verifyEqual(redSys.dim, 10);
        end

        function MIMOWithFixedTangentDir(testCase,sysMIMO)
            % Test MIMO system reduction with initial shifts supplied
            redSys = trksmPH(sysMIMO,8,1e-3,1e2,2);
            testCase.verifyEqual(redSys.dim, 8);

            % Test MIMO system reduction with initial shifts supplied
            redSys = trksmPH(sysMIMO,4,1e-3,1e2,1);
            testCase.verifyEqual(redSys.dim, 4);

        end

    end

    methods (Test, TestTags = {'Output'})
        % check (optional) output

        function interpolationConditions(testCase, sys, sysMIMO)
            % Verify that the ROM interpolates the FOM at the returned shift vector and tangent
            % directions
            % SISO
            [redSys] = trksmPH(sys,8);
            for i=1:length(redSys.info.shifts)
                testCase.verifyLessThan(norm(sys.evalfr(redSys.info.shifts(i))-redSys.evalfr(redSys.info.shifts(i))),testCase.tolerance)
            end
            % MIMO
            [redSys] = trksmPH(sysMIMO,8);
            for i=1:length(redSys.info.shifts)
                testCase.verifyLessThan(norm(sysMIMO.evalfr(redSys.info.shifts(i))*redSys.info.Dr{i}-redSys.evalfr(redSys.info.shifts(i))*redSys.info.Dr{i}),testCase.tolerance)
            end
        end

        function redSysInfo(testCase, sys)
            % Verify that redSys object has additional information
            redSys = trksmPH(sys,4);

            testCase.verifyTrue(isstruct(redSys.info));
            testCase.verifyTrue(isfield(redSys.info, "residual"));
            testCase.verifyTrue(isfield(redSys.info, "shifts"));
            testCase.verifyTrue(isfield(redSys.info, "Dr"));
        end
    end

    methods (Test, TestTags = {'Options'})
        % check if options are applied correctly

        function options_isReal(testCase,sysMIMO)
            % Check if real shifts are chosen depending on opts.isReal
            opts.isReal = true;
            redSys = trksmPH(sysMIMO,6,opts);
            testCase.verifyTrue(isreal(redSys.info.shifts));
            testCase.verifyEqual(redSys.parameters.isReal,true);

            opts.isReal = false;
            redSys = trksmPH(sysMIMO,6,opts);
            testCase.verifyTrue(~isreal(redSys.info.shifts));
            testCase.verifyTrue(~redSys.parameters.isReal);
        end

        function options_border(testCase,sys)
            % Check option opts.border
            opts.border = true;
            redSys = trksmPH(sys,2,opts);
            testCase.verifyTrue(redSys.parameters.border);
            opts.border = false;
            redSys = trksmPH(sys,2,opts);
            testCase.verifyTrue(~redSys.parameters.border);
        end

        function options_log(testCase,sys)
            % Check option opts.log
            opts.log = true;
            redSys = trksmPH(sys,2,opts);
            testCase.verifyTrue(redSys.parameters.log);
            opts.log = false;
            redSys = trksmPH(sys,2,opts);
            testCase.verifyTrue(~redSys.parameters.log);
        end

        function options_tol(testCase,sys)
            % Check option opts.tol
            opts.tol = Inf;
            redSys = trksmPH(sys,10,opts);
            testCase.verifyLessThan(redSys.dim,10);
        end

    end

    methods (Test, TestTags = {'Functionality'})
        % Verify correct outputs for more realistic scenarios

        function accuracy(testCase)
            % Check accuracy for order 200 model

            sysLarge = setup_MassSpringDamperSystem(200,2,1,1);

            Opts.tol = 1e-8;
            redSys = trksmPH(sysLarge,100,Opts);
            testCase.verifyLessThan(norm(redSys-sysLarge)/norm(sysLarge), 1e-5);
        end
    end

end