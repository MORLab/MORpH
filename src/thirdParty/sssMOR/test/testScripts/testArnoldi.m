classdef testArnoldi < sssTest
% testArnoldi - testing of arnoldi.m
%
% Description:
%   The function arnoldi.m is tested (5 tests) on:
%    + comparing the columns of V to a directly calculated solution 
%    + AV=A*V, EV=E*V
%    + V'*V orthonormal
%    + V purely real
%    + rank(V) full
%    + Neither Inf nor NaN in V
%    + s0: purely real, purely imaginary, zero, Inf (Markov-parameter)
%    + test systems: building, iss, fom, eady, rail_1357 (with E-matrix)
%
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
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    
    methods(Test)        
        function testArnoldi1(testCase)
            %s0=Inf, real, imag (without E-matrix)
            load('building.mat');

            [V] = arnoldi(speye(size(A)),A,B,[Inf, 50, 100, 200, 300, 1-1i, 1+1i]);
            actSolution={V};
            
            temp=full(real(B));
            expV1=temp/norm(temp);      
            
            temp=(A-50*eye(size(A)))\B;
            expV2=temp/norm(temp);
            
            temp=(A-100*eye(size(A)))\B;
            expV3=temp/norm(temp);
            
            temp=(A-200*eye(size(A)))\B;
            expV4=temp/norm(temp);
            
            temp=(A-300*eye(size(A)))\B;
            expV5=temp/norm(temp);
            
            temp=real((A-(1-1i)*eye(size(A)))\B);
            expV6=temp/norm(temp);
            
            temp=imag((A-(1-1i)*eye(size(A)))\B);
            expV7=temp/norm(temp);
     
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5, expV6, expV7]);
            expV=[-expV(:,1),expV(:,2), expV(:,3), -expV(:,4), -expV(:,5), -expV(:,6), expV(:,7)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testArnoldi2(testCase) 
            %s0: only imag (multiple value)
            sys = sss('iss.mat');
            sys = sys(1,1);

            [V] = arnoldi(speye(size(sys.A)),sys.A,sys.B,[1-1i, 1-1i, 1-1i, 1+1i, 1+1i, 1+1i]);
            actSolution={V};
            
            temp=real((sys.A-(1-1i)*eye(size(sys.A)))\sys.B);
            expV1=temp/norm(temp);
 
            temp=real((sys.A-(1-1i)*eye(size(sys.A)))\expV1);
            expV2=temp/norm(temp);
            
            temp=real((sys.A-(1-1i)*eye(size(sys.A)))\expV2);
            expV3=temp/norm(temp);
            
            temp=imag((sys.A-(1-1i)*eye(size(sys.A)))\sys.B);
            expV4=temp/norm(temp);
            
            temp=imag((sys.A-(1-1i)*eye(size(sys.A)))\expV1);
            expV5=temp/norm(temp);
            
            temp=imag((sys.A-(1-1i)*eye(size(sys.A)))\expV2);
            expV6=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5, expV6]);
            expV=[expV(:,1), -expV(:,2), expV(:,3), -expV(:,4),expV(:,5), -expV(:,6)];
            
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        function testArnoldi3(testCase) 
            %s0: only real (multiple value)
            load('eady.mat');

            [V] = arnoldi(speye(size(A)),A,B,[50,50,50,50]);
            actSolution={V};
     
            temp=(A-50*eye(size(A)))\B;
            expV1=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV1;
            expV2=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV2;
            expV3=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV3;
            expV4=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4]);
            expV=[-expV(:,1), -expV(:,2), expV(:,3), expV(:,4)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testArnoldi4(testCase) 
            %s0: only Inf (multiple value)
            load('fom.mat');

            [V] = arnoldi(speye(size(A)),A,B,[Inf, Inf, Inf, Inf]);
            actSolution={V};
     
            temp=full(B);
            expV1=temp/norm(temp);
            
            temp=full(A*B);
            expV2=temp/norm(temp);

            temp=full(A^2*B);
            expV3=temp/norm(temp);
            
            temp=full(A^3*B);
            expV4=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4]);
            expV=[-expV(:,1), -expV(:,2), -expV(:,3), expV(:,4)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        function testArnoldi5(testCase) 
            %s0=Inf, real, imag, zero with E-matrix
            sys = sss('rail_1357.mat');
            sys = sys(1,1); E = sys.E; A = sys.A; B = sys.B;
            
            [V] = arnoldi(E,A,B,[Inf, 0, 100, 4+13i, 4-13i], @(x,y) (x'*y));
            actSolution={V};
             
            temp=full(real(E\B));
            expV1=temp/norm(temp);      
            
            temp=full((A-0*E)\B);
            expV2=temp/norm(temp);
            
            temp=full((A-100*E)\B);
            expV3=temp/norm(temp);
            
            temp=full(real((A-(4+13i)*E)\B));
            expV4=temp/norm(temp);
            
            temp=full(imag((A-(4+13i)*E)\B));
            expV5=temp/norm(temp);
     
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5]);
            expV=[-expV(:,1), expV(:,2), expV(:,3), expV(:,4), expV(:,5)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testArnoldi6(testCase) 
            %test Hermite arnoldi for SISO systems
            for i=1:length(testCase.sysCell)
                %  test system
                sys=testCase.sysCell{i};
                sys = sys(1,1);
               
                %  get good shifts
                n = 6; r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r, l);
                Opts.rType = 'dir';
                [r,p] = residue (sysrIrka,Opts);
                s0 = -(conj(p)).'; Lt = r{1}; Rt = r{2}.';         
                
                %   run Hermite arnoldi
                [V, ~, Rv, W, ~, Lw] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0, Rt, Lt,@(x,y) (x'*y));
                actSolution={W};
                
                %   run output arnoldi
                [Wexp,~, LwExp] = arnoldi(sys.E.',sys.A.',sys.C.',s0, Lt, @(x,y) (x'*y));
                expSolution= {Wexp};
                
                %   Verify W
                verification(testCase, actSolution,expSolution,W)
                
                %   Verify Lsylv equality
                verifyEqual(testCase, Lw, LwExp, 'RelTol', 1e-7,...
                    'Generalized tangential directions do not match');
                
               %   Verify solution of Sylvester EQ
               %       AV - EV(Er\Ar) - B_R = 0
               %       A.'W - E.'W (Er.'\Ar.') - C_ L = 0
               sysr = sss(W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V);
               B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
               res1 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rv);
               
               % Rexp = (B_.'*B_)\B_.'*(sys.A*V - sys.E*V*(sysr.E\sysr.A));
               % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rexp)
               
               sysd = sys.'; sysrd = sysr.'; %dual systems
               C_ = sysd.B - sysd.E*W*(sysrd.E\sysrd.B);
               res2 = norm(sysd.A*W - sysd.E*W*(sysrd.E\sysrd.A) - C_*Lw);
               
               verifyEqual(testCase, [res1, res2] , [0, 0], 'AbsTol', 1e-7,...
                    'Sylvester EQ is not satisfied');
                
               % verify moment matching
              actM = moments(sysr,s0, 2); expM = moments(sys,s0, 2);
              actMt = {}; expMt = {};
              for iS = 1:length(s0)
                  actMt = [actMt, {actM(:,:,iS)*Rt(:,iS), Lt(:,iS).'*actM(:,:,iS),...
                           Lt(:,iS).'*actM(:,:,iS)*Rt(:,iS), Lt(:,iS).'*actM(:,:,iS+1)*Rt(:,iS)}];
                  expMt = [expMt, {expM(:,:,iS)*Rt(:,iS), Lt(:,iS).'*expM(:,:,iS),...
                           Lt(:,iS).'*expM(:,:,iS)*Rt(:,iS), Lt(:,iS).'*expM(:,:,iS+1)*Rt(:,iS)}];
              end
              verifyEqual(testCase, actMt, expMt, 'RelTol', 1e-6);
            end
        end
        function testArnoldi7(testCase) 
            %test Hermite arnoldi for MIMO systems
            for i=1:length(testCase.sysCell)
                %  test system
                sys=testCase.sysCell{i};
                if ~sys.isSiso
                %  get good shifts
                n = 6; r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r, l);
                Opts.rType = 'dir';
                [r,p] = residue (sysrIrka,Opts);
                s0 = -(conj(p)).'; Lt = r{1}; Rt = r{2}.';         
                
                %   run Hermite arnoldi
                [V,~,Rv,W,~,Lw] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0, Rt, Lt,@(x,y) (x'*y));
                actSolution={W}; sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V,sys.D,W'*sys.E*V);
                
                %   run output arnoldi
                [Wexp,~, LwExp] = arnoldi(sys.E.',sys.A.',sys.C.',s0, Lt, @(x,y) (x'*y));
                expSolution= {Wexp};
                
                %   Verify W
                verification(testCase, actSolution,expSolution,W)
                
                %   Verify Lsylv equality vs onsided
                verifyEqual(testCase, Lw, LwExp, 'RelTol', 1e-7,...
                    'Generalized tangential directions do not match');
                %   Verify Rsylv, Lsylv vs getSylvester
                [expRsylv] = getSylvester(sys,sysr,V);
                [expLsylv] = getSylvester(sys,sysr,W,'W');
              
                verifyEqual(testCase, {Rv, Lw}, {expRsylv,expLsylv}, ...
                    'RelTol', 1e-7, 'Generalized tangential directions do not match');
                
               %   Verify solution of Sylvester EQ
               %       AV - EV(Er\Ar) - B_R = 0
               %       A.'W - E.'W (Er.'\Ar.') - C_ L = 0
               sysr = sss(W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V);
               B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
               res1 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rv);
               
               % Rexp = (B_.'*B_)\B_.'*(sys.A*V - sys.E*V*(sysr.E\sysr.A));
               % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rexp)
               
               sysd = sys.'; sysrd = sysr.'; %dual systems
               C_ = sysd.B - sysd.E*W*(sysrd.E\sysrd.B);
               res2 = norm(sysd.A*W - sysd.E*W*(sysrd.E\sysrd.A) - C_*Lw);
               
               verifyEqual(testCase, [res1, res2] , [0, 0], 'AbsTol', 1e-7,...
                    'Sylvester EQ is not satisfied');
                
              % verify moment matching
              actM = moments(sysr,s0, 2); expM = moments(sys,s0, 2);
              actMt = {}; expMt = {};
              for iS = 1:length(s0)
                  jM = iS*2-1;
                  actMt = [actMt, {actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM),...
                           Lt(:,iS).'*actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM+1)*Rt(:,iS)}];
                  expMt = [expMt, {expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM),...
                           Lt(:,iS).'*expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM+1)*Rt(:,iS)}];
              end
              verifyEqual(testCase, actMt, expMt, 'RelTol', 1e-6);
              end
            end
        end
        function testArnoldiMIMORail(testCase) 
            %test Hermite arnoldi for rail system
            sys = sss('rail_1357');
            sysCell{1} = sys;
            sysCell{2} = sys(1:sys.p,1);
            sysCell{3} = sys(1,1:sys.m);
            sysCell{4} = sys(1:sys.p,1:sys.p);
            n = 6;
            for iCase=1:4
                sys = sysCell{iCase};
                r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r, l);
                Opts.rType = 'dir';
                [r,p] = residue (sysrIrka,Opts);
                s0 = -(conj(p)).'; Lt = r{1}; Rt = r{2}.';         

                %   run Hermite arnoldi
                [V,~,Rv,W,~,Lw] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0, Rt, Lt,@(x,y) (x'*y));
                actSolution={W}; sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V,sys.D,W'*sys.E*V);

                %   run output arnoldi
                [Wexp,~, LwExp] = arnoldi(sys.E.',sys.A.',sys.C.',s0, Lt, @(x,y) (x'*y));
                expSolution= {Wexp};

                %   Verify W
                verification(testCase, actSolution,expSolution,W)

                %   Verify Lsylv equality vs onsided
                verifyEqual(testCase, Lw, LwExp, 'RelTol', 1e-7,...
                    'Generalized tangential directions do not match');
                %   Verify Rsylv, Lsylv vs getSylvester
                [expRsylv] = getSylvester(sys,sysr,V);
                [expLsylv] = getSylvester(sys,sysr,W,'W');

                verifyEqual(testCase, {Rv, Lw}, {expRsylv,expLsylv}, ...
                    'RelTol', 1e-7, 'Generalized tangential directions do not match');

               %   Verify solution of Sylvester EQ
               %       AV - EV(Er\Ar) - B_R = 0
               %       A.'W - E.'W (Er.'\Ar.') - C_ L = 0
               sysr = sss(W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V);
               B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
               res1 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rv);

               % Rexp = (B_.'*B_)\B_.'*(sys.A*V - sys.E*V*(sysr.E\sysr.A));
               % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rexp)

               sysd = sys.'; sysrd = sysr.'; %dual systems
               C_ = sysd.B - sysd.E*W*(sysrd.E\sysrd.B);
               res2 = norm(sysd.A*W - sysd.E*W*(sysrd.E\sysrd.A) - C_*Lw);

               verifyEqual(testCase, [res1, res2] , [0, 0], 'AbsTol', 1e-7,...
                    'Sylvester EQ is not satisfied');

              % verify moment matching
              actM = moments(sysr,s0, 2); expM = moments(sys,s0, 2);
              actMt = {}; expMt = {};
              for iS = 1:length(s0)
                  jM = iS*2-1;
                  actMt = [actMt, {actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM),...
                           Lt(:,iS).'*actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM+1)*Rt(:,iS)}];
                  expMt = [expMt, {expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM),...
                           Lt(:,iS).'*expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM+1)*Rt(:,iS)}];
              end
              verifyEqual(testCase, actMt, expMt, 'AbsTol', 1e-6);
            end
        end
        function testArnoldi8(testCase) 
            %test Sylvester matrices
            Opts.orth='2mgs';
            Opts.reorth='qr';
            Opts.real='real';
            s0=[50, 100, 200, 300, 1-1i, 1+1i];
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                [V, Sv, Rv] = arnoldi(sys.E,sys.A,sys.B,s0,Opts);
                % residual of sylvester equation
                actSolution=norm(sys.A*V-sys.E*V*Sv-sys.B*Rv);
                verifyLessThan(testCase,actSolution,1e-5);
            end
        end
    end
end

function [] = verification(testCase, actSolution,expSolution,V)
%        verifyEqual(testCase, actSolution, expSolution,'RelTol',0.2,'AbsTol',0.000000001,...
%             'Difference between actual and expected exceeds relative tolerance');
       s = svd([actSolution{:},expSolution{:}]); s = s/s(1); 
       verifyEqual(testCase, {sum(s>1e-6)}, {size(actSolution{:},2)},...
           'Arnoldi failed to recreate the correct subspace')
       verifyLessThanOrEqual(testCase, norm(V'*V-speye(size(V,2)),'fro'), ...
           1e-12, 'Orthonormalization failed');
       verifyLessThanOrEqual(testCase, isreal(V), 1, ...
            'V is not purely real'); 
       verifyEqual(testCase, rank(V), size(V,2),...
            'Rank(V) is not full');
       verifyEqual(testCase, nnz(isinf(V)), 0, ...
            'V contains Inf');
       verifyEqual(testCase, nnz(isnan(V)), 0, ...
            'V contains Nan');
        verifyEqual(testCase, size(actSolution), size(expSolution),...
            'V is not of the correct size')
        
end