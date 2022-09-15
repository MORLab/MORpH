classdef testIsH2opt < sssTest
% testIsH2opt - testing of isH2opt.m
%
% Description:
%   The function isH2opt.m is tested (2 tests) on:
%    + is H2opt: s0 mirrored and moments matching
%    + not H2opt
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  28 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
    
    methods(Test)
        function testIsH2opt1(testCase)
            %local H2-optimal
            load('building.mat');
            opts=struct('tol',0.0005,'maxiter',300,'type','stab','stopCrit','s0');

            s0=[0,0,50,100,Inf, Inf, 1+5i,1-5i,14-0.2i,14+0.2i];
            sys=sss(A,B,C);
            [sysr, ~, ~, s0] = irka(sys,s0,opts);
             
            opts=struct('tol',0.0005);
            opt = isH2opt(sys,sysr,s0,opts);
            actSolution={opt};
            
            %test if s0 is mirrored
            verifyEqual(testCase, sort(s0.'), sort(-(conj(eig(sysr)))),...
                'RelTol',0.1,'AbsTol',0.000001,'s0 not mirrored');
            
            %test if moments are matching
            for i=1:length(s0)
            verifyEqual(testCase, moments(sys,s0(i),2),  moments(sysr,s0(i),2),...
                'RelTol',0.1,'AbsTol',0.000001,'moments not matching');
            end

            expSolution={true};
            verification(testCase, actSolution, expSolution);
        end
        function testIsH2opt2(testCase)
            %not local H2-optimal ('Warning: IRKA has not converged after 3
            %steps')
            load('building.mat');
            opts=struct('epsilon',0.1,'maxiter',3,'type','stab','stopCrit','sysr');
            
            s0=[0,0,50,100,Inf, Inf, 1+5i,1-5i,14-0.2i,14+0.2i];
            sys=sss(A,B,C);
            [sysr, ~, ~, s0] = irka(sys,s0,opts);
             
            opts=struct('tol',0.005);
            opt = isH2opt(sys,sysr,s0,opts);
            actSolution={opt};
            
            expSolution={false};
            verification(testCase, actSolution, expSolution);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,...
            'Difference between actual and expected exceeds relative tolerance');
end