classdef test_testIntro < matlab.unittest.TestCase
% Tests must inherit from class matlab.unittest.TestCase
% You can run a testscript like this one by pressing <F5> on your keyboard,
% selecting "Run Tests" in the navigation bar (see "Editor"), or by
% executing "runtests("test_testIntro")" in the command window.
% During testing, test-triggered outputs will be shown in the command
% window. In addition, you will be given a "TestResult"-Object in the end
% which contains a summary of all tests on whether they passed or failed.

%% As for any class, you can specify properties of objects:
    
    properties (TestParameter)
    % Specify the parameters for your test here. Each parameter must be a
    % cell array containing all parameters to test.
        param1 = {1,2,3};
        param2 = {"1","2","3"};
    end
    
    properties
    % Specify reused variables which are not explicit test parameters here.
        tolerance = 1;
    end
    
%% If you want to set up a "clean" testing environment before running tests, use fixture functions
%   see also: https://www.mathworks.com/help/matlab/matlab_prog/write-setup-and-teardown-code-using-classes.html
    
    methods (TestMethodSetup)
        function lessOutput(testCase)
        % example: disable all warnings
            warning('off','all');
            testCase.addTeardown(@warning,'on','all');
        end
    end
    
    methods (TestMethodTeardown)
        function activateWarnings(testCase)
            % warning('on','all');
            % better: Include testCase.addTeardown(@warning,'on','all');
        end
    end
    
%% Write your unit tests within special Test-method-blocks:
% Any error within a test function will stop the execution of this
% particular function but not the test as a whole. Testing will simply
% continue with the next function.

    methods (Test)
        function testfunction1(testCase)
        % test functions always need a testCase-object as first input (name
        % is arbitrary). Test parameters may be used optionally.
            
            % You can perform any computations within test functions
            x = 1;
            x = x+1;
            %etc...
            
            % To verify results, you could use the assert function
            assert(x == 2);
            % which will throw an error if the statement is incorrect.
            
            % However, we recommend to use the
            % matlab.unittest.qualifications package instead, see:
            % "https://www.mathworks.com/help/matlab/matlab_prog/types-of-qualifications.html"
            % example:
            testCase.verifyEqual(x,2)       % will not prompt failed verification
            testCase.verifyLessThan(x,1)    % will prompt failed verification
            
            % These function will not stop the test function if the
            % statement is incorrect (as "assert" would). In addition, you
            % will get more detailed information on the command window.
        end
    end
    
    % Use TestTags = {...} to group tests together. You can use this
    % later to only run tests with a certain tag. You can also assign
    % multiple test tags to one block of methods by extending the
    % respective cell array in the method header.
    methods (Test, TestTags = {'Tag1'})
        % Be including test parameters as input parameters, you may use the
        % cell entries within any test function.
        function testfunction2_a(testCase,param1)
            testCase.verifyEqual(param1,2);
        end
        % You can also include other class properties by using the testCase-object:
        function testfunction2_b(testCase, param1)
            testCase.verifyLessThanOrEqual(abs(param1-1),testCase.tolerance)
        end
    end
    
    % If you have multiple test parameters, you can specify how Matlab
    % combines them:
    % Use "ParameterCombination = 'Sequential'" to iterate through
    % parameter combinations with same index.
    methods (Test, TestTags = {'Tag2'},ParameterCombination = 'sequential')
        function testfunction3(testCase,param1,param2)
            testCase.verifyEqual(param1,str2double(param2));
        end
    end
    % Use "ParameterCombination = 'Exhaustive'" to iterate through all
    % possible combinations of parameters.
    % See "https://www.mathworks.com/help/matlab/matlab_prog/create-advanced-parameterized-test.html"
    % for more parameter combinations.
    methods (Test, TestTags = {'Tag2'},ParameterCombination = 'exhaustive')
        function testfunction4(testCase,param1,param2)
            testCase.verifyLessThanOrEqual(param1,str2double(param2));
        end
    end

end