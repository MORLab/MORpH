classdef test_testTemplate < matlab.unittest.TestCase
% TEST testName - Tests the functionalities of <module name>
%
% Press <F5> or enter "runtests("test_testName")" to run this testscript
%
% Description:
%   < Describe what this test is intended to do here. Include, e.g.,
%     considered test scenarios (grouped by tags)                    >
%
%-----------------------------------------------------------------------
% This file is part of 
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      <Name of authors>
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a> 
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

%% Properties

    properties (TestParameter) % test parameters
        
    end
    
    properties % other properties (not used as test inputs)
        
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
    
        function basicInput(testCase)
        % Describe what this test is supposed to verify, e.g.
        % Verify that 1 equals to 2
        
            testCase.verifyEqual(1,2)
        end
        
    end
    
    methods (Test, TestTags = {'Output'})
    % Check if (optional) output is correct
    
    end
    
    methods (Test, TestTags = {'Options'})
    % Check if options are applied correctly
    
    end

    methods (Test, TestTags = {'Functionality'})
    % Check if functionality is given (e.g. accuracy of results)
    
    end
    
    methods (Test, TestTags = {'Benchmarks'})
    % Use benchmarks if available
    
    end
    
%% Supporting functions (e.g. for reuse)

    methods
        
    end

end