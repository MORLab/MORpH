function [dPhi_dx,dLambda_dx] = dPR_dx(sys,theta,dims,evp,method)
% dPR_dx - Derivative of poles Lambda and residues Phi with respect to
% port-Hamiltonian matrix entries
%
% Syntax:
%   [dPhi_dx,dLambda_dx] = dPR_dx(sys,theta,dims,evp)
%
% Description:
%       Computes the derivative of poles and residues with respect to
%       port-Hamiltonian matrices (E,J,R,Q,G,P) and Lyapunov factor L%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object of the model
%       - theta:    current parameter vector (optimization point)
%       - dims:     Dimensions of subvectors in theta
%       - evp:      Struct with matrices of eigenvalue problem
%                       A*evp.Xr = evp.Xr*evp.Lambda
%                       evp.Xl'*A = evp.Lambda*evp.Xl'
%       - method:   Method for differentiating the eigenvalue problem
%                   ['magnus','murthy','rudisill']
%
% Output Arguments:
%       - dPhi_dx:      Derivative of residual w.r.t. to theta
%       - dLambda_dx: 	Derivative of poles w.r.t. to theta
%
%
% See Also:
%       prOpt
%
% References:
%       [1] J. R. Magnus, "On differentiating eigenvalues and
%           eigenvectors”, Econometric Theory, vol. 1, no. 02, pp. 179–191, 1985.
%       [2] D. V. Murthy and R. T. Haftka, “Derivatives of eigenvalues and
%           eigenvectors of a general complex matrix”, International Journal for
%           Numerical Methods in Engineering, vol. 26, no. 2, pp. 293–311, 1988.
%       [3] C. S. Rudisill, “Derivatives of eigenvalues and eigenvectors
%           for a general matrix" , AIAA Journal, vol. 12, no. 5, pp. 721–722, 1974.
%       [4] P. Schwerdtner et al. Structure-Preserving Model Order
%           Reduction for Index One Port-Hamiltonian Descriptor Systems. 
%           arXiv Preprint arXiv:2206.01608. 2022. url: https://arxiv.org/abs/2206.01608.
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

%% Initialization
% Dimensions
r = dims.nx;
m = dims.nu;

% Eigenvalue problem
lambda = diag(evp.Lambda);
Xr = evp.Xr;
Xl = evp.Xl;

% State-space matrices
E = sys.E;
A = (sys.J-sys.R)*sys.Q;
B = sys.G-sys.P;
C = (sys.G+sys.P)'*sys.Q;

% Factorization of sys.S
Ls = zeros(m,m);
if ~isequal(sys.S,zeros(m,m))
    [Us,Ds] = eig(full(sys.S));
    Ls = Us*sqrt(Ds);
end

% Some frequently used matrices
Erinv = E\eye(r);
Xrinv = Xr\eye(r);
EXinv = (Xrinv*Erinv).';
EXinvB = Xrinv*Erinv*B;
indMat = reshape(1:r*r,[r,r]);

% Scale EVP in case of 'murthy'
switch method
    case 'murthy'
        % Determine row m of max element
        [~,m_max] = max(abs(Xrinv).'.*abs(Xr),[],1);
        Gamma = diag(1./diag(Xr(m_max,:)));
        % Transformation of EVP
        Xr = Xr*Gamma;
    otherwise
end

% Pre-allocate matrices
if dims.nE,[dLambda_dE,dPhi_dE,dl_dE,dr_dE,dXr_dE_rsh] = prealloc(r,m,dims.nE); end
if dims.nJ,[dLambda_dJ,dPhi_dJ,dl_dJ,dr_dJ,dXr_dJ_rsh] = prealloc(r,m,dims.nJ); end
if dims.nR,[dLambda_dR,dPhi_dR,dl_dR,dr_dR,dXr_dR_rsh] = prealloc(r,m,dims.nR); end
if dims.nQ,[dLambda_dQ,dPhi_dQ,dl_dQ,dr_dQ,dXr_dQ_rsh] = prealloc(r,m,dims.nQ); end
if dims.nP,[dLambda_dP,dPhi_dP,dl_dP,dr_dP,dXr_dP_rsh] = prealloc(r,m,dims.nP); end
if dims.nG,[dLambda_dG,dPhi_dG] = prealloc(r,m,dims.nG); end
if dims.nL,[dLambda_dL,dPhi_dL] = prealloc(r,m,dims.nL); end

% Get sub-vectors of theta
thetaE = theta(1:dims.nE); theta = theta(dims.nE+1:end);
theta = theta(dims.nJ+1:end);
thetaR = theta(1:dims.nR); theta = theta(dims.nR+1:end);
thetaQ = theta(1:dims.nQ); theta = theta(dims.nQ+1:end);
theta = theta(dims.nG+1:end);
thetaP = theta(1:dims.nP);
thetaL = theta(1:dims.nL);

%% Derivatives w.r.t. matrices where dLambdaR_dx ~= 0
ptr = 1;
for j=1:r
    switch method
        case 'rudisill'
            C_rud = [Xr(:,j)',0;A-lambda(j)*E,-E*Xr(:,j)];
        case 'murthy'
            C_rud = [A-lambda(j)*E,-E*Xr(:,j)];
            C_rud = C_rud(:,[1:m_max(j)-1,m_max(j)+1:end]);
        otherwise
    end

    % E
    if dims.nE
        % Block formulation
        switch method
            case 'rudisill'
                D_rud = [zeros(1,r*r);lambda(j)*(kron(vtu(thetaE)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaE)',eye(r)))];
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                dLambda_dE(j,:) = slv(end,:);
                dZj_dE = slv(1:end-1,:);
                dl_dE(ptr:ptr+m-1,:) = C*dZj_dE;
                dXr_dE_rsh(:,j) = dZj_dE(:);
            case 'murthy'
                D_rud = lambda(j)*(kron(vtu(thetaE)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaE)',eye(r)));
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                slv = [slv(1:m_max(j)-1,:);zeros(1,size(slv,2));slv(m_max(j):end,:)];
                dLambda_dE(j,:) = slv(end,:);
                dZj_dE = slv(1:end-1,:);
                dl_dE(ptr:ptr+m-1,:) = C*dZj_dE;
                dXr_dE_rsh(:,j) = dZj_dE(:);
            case 'magnus'
                dZj_dE = -pinv(lambda(j)*E-A,1e-12)*(lambda(j)*eye(r,r)-(E*Xr(:,j)*Xl(:,j)'*lambda(j))/(Xl(:,j)'*E*Xr(:,j)))*(kron(vtu(thetaE)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaE)',eye(r)));
                % Remove parameters which are zero
                ind = nonzeros(tril(indMat));
                dZj_dE = dZj_dE(:,ind);
                dl_dE(ptr:ptr+m-1,:) = C*dZj_dE;
                dXr_dE_rsh(:,j) = dZj_dE(:);

                dLambdaRi_dE = (-Xl(:,j)'*lambda(j)*(kron(vtu(thetaE)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaE)',eye(r))))/(Xl(:,j)'*E*Xr(:,j));
                dLambdaRi_dE = dLambdaRi_dE(:,ind);
                dLambda_dE(j,:) = dLambdaRi_dE(:).';

            otherwise
        end
    end

    % J
    if dims.nJ
        % Block formulation
        switch method
            case 'rudisill'
                D_rud = [zeros(1,r*r);-kron((sys.Q*Xr(:,j)).',eye(r))];
                % Remove parameters which are zero
                ind1 = tril(indMat,-1);ind1 = nonzeros(ind1);
                ind2 = triu(indMat,1)';ind2 = nonzeros(ind2);
                D_rud = D_rud(:,ind1)-D_rud(:,ind2);
                slv = C_rud\D_rud;
                dLambda_dJ(j,:) = slv(end,:);
                dZj_dJ2 = slv(1:end-1,:);
                dl_dJ(ptr:ptr+m-1,:) = C*dZj_dJ2;
                dXr_dJ_rsh(:,j) = dZj_dJ2(:);
            case 'murthy'
                D_rud = -kron((sys.Q*Xr(:,j)).',eye(r));
                % Remove parameters which are zero
                ind1 = tril(indMat,-1);ind1 = nonzeros(ind1);
                ind2 = triu(indMat,1)';ind2 = nonzeros(ind2);
                D_rud = D_rud(:,ind1)-D_rud(:,ind2);
                slv = C_rud\D_rud;
                % Set derivative of fixed eigenvector entries to 0
                slv = [slv(1:m_max(j)-1,:);zeros(1,size(slv,2));slv(m_max(j):end,:)];
                dLambda_dJ(j,:) = slv(end,:);
                dZj_dJ2 = slv(1:end-1,:);
                dl_dJ(ptr:ptr+m-1,:) = C*dZj_dJ2;
                dXr_dJ_rsh(:,j) = dZj_dJ2(:);
            case 'magnus'
                dZj_dJ = pinv(lambda(j)*E-A,1e-12)*(eye(r,r)-(E*Xr(:,j)*Xl(:,j)')/(Xl(:,j)'*E*Xr(:,j)))*kron((sys.Q*Xr(:,j)).',eye(r,r));
                % Remove parameters which are zero
                ind1 = tril(indMat,-1);ind1 = nonzeros(ind1);
                ind2 = triu(indMat,1)';ind2 = nonzeros(ind2);
                dZj_dJ = dZj_dJ(:,ind1)-dZj_dJ(:,ind2);
                % Construct reshaped version of dZ_dJr for dR_dJr
                dXr_dJ_rsh(:,j) = dZj_dJ(:);
                dl_dJ(ptr:ptr+m-1,:) = C*dZj_dJ;
                dLambdaRi_dJ = (sys.Q*Xr(:,j)*Xl(:,j)').'/(Xl(:,j)'*E*Xr(:,j));
                dLambdaRi_dJ = dLambdaRi_dJ(:).';
                dLambdaRi_dJ = dLambdaRi_dJ(:,ind1)-dLambdaRi_dJ(:,ind2);
                dLambda_dJ(j,:) = dLambdaRi_dJ(:).';
            otherwise
        end
    end

    % R
    if dims.nR
        % Block formulation
        switch method
            case 'rudisill'
                D_rud = [zeros(1,r*r);kron(vtu(thetaR)',(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtu(thetaR)',eye(r))];
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                dLambda_dR(j,:) = slv(end,:);
                dZj_dR = slv(1:end-1,:);
                dl_dR(ptr:ptr+m-1,:) = C*dZj_dR;
                dXr_dR_rsh(:,j) = dZj_dR(:);
            case 'murthy'
                D_rud = kron(vtu(thetaR)',(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtu(thetaR)',eye(r));
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                slv = [slv(1:m_max(j)-1,:);zeros(1,size(slv,2));slv(m_max(j):end,:)];
                dLambda_dR(j,:) = slv(end,:);
                dZj_dR = slv(1:end-1,:);
                dl_dR(ptr:ptr+m-1,:) = C*dZj_dR;
                dXr_dR_rsh(:,j) = dZj_dR(:);
            case 'magnus'
                dZj_dR = -pinv(lambda(j)*E-A,1e-12)*(eye(r,r)-(E*Xr(:,j)*Xl(:,j)')/(Xl(:,j)'*E*Xr(:,j)))*(kron(vtu(thetaR)',(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtu(thetaR)',eye(r)));
                % Remove parameters which are zero
                ind = nonzeros(tril(indMat));
                dZj_dR = dZj_dR(:,ind);
                dl_dR(ptr:ptr+m-1,:) = C*dZj_dR;
                dXr_dR_rsh(:,j) = dZj_dR(:);

                dLambdaRi_dR = -(Xl(:,j)'*(kron(vtu(thetaR)',(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtu(thetaR)',eye(r))))/(Xl(:,j)'*E*Xr(:,j));
                dLambdaRi_dR = dLambdaRi_dR(:,ind);
                dLambda_dR(j,:) = dLambdaRi_dR(:).';
            otherwise
        end
    end

    % Q
    if dims.nQ
        % Block formulation
        switch method
            case 'rudisill'
                D_rud = [zeros(1,r*r);-(sys.J-sys.R)*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)))];
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                dLambda_dQ(j,:) = slv(end,:);
                dZj_dQ = slv(1:end-1,:);

                % Block formulation
                slv2 = (sys.G+sys.P)'*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)));
                dl_dQ(ptr:ptr+m-1,:) = C*dZj_dQ + slv2(:,ind);

                dXr_dQ_rsh(:,j) = dZj_dQ(:);
            case 'murthy'
                D_rud = -(sys.J-sys.R)*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)));
                ind = nonzeros(tril(indMat));
                D_rud = D_rud(:,ind);
                slv = C_rud\D_rud;
                % Set derivative of fixed eigenvector entries to 0
                slv = [slv(1:m_max(j)-1,:);zeros(1,size(slv,2));slv(m_max(j):end,:)];
                dLambda_dQ(j,:) = slv(end,:);
                dZj_dQ = slv(1:end-1,:);

                % Block formulation
                slv2 = (sys.G+sys.P)'*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)));
                dl_dQ(ptr:ptr+m-1,:) = C*dZj_dQ + slv2(:,ind);

                dXr_dQ_rsh(:,j) = dZj_dQ(:);

            case 'magnus'
                dZj_dQ = pinv(lambda(j)*E-A,1e-12)*(eye(r,r)-(E*Xr(:,j)*Xl(:,j)')/(Xl(:,j)'*E*Xr(:,j)))*(sys.J-sys.R)*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)));
                % Remove parameters which are zero
                ind = nonzeros(tril(indMat));
                dZj_dQ = dZj_dQ(:,ind);
                % Block formulation
                slv2 = (sys.G+sys.P)'*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r)));
                dl_dQ(ptr:ptr+m-1,:) = C*dZj_dQ + slv2(:,ind);
                dXr_dQ_rsh(:,j) = dZj_dQ(:);

                dLambdaRi_dQ = (Xl(:,j)'*(sys.J-sys.R)*(kron(vtu(thetaQ)',(Xr(:,j)).')+kron((Xr(:,j)).'*vtu(thetaQ)',eye(r))))/(Xl(:,j)'*E*Xr(:,j));
                dLambdaRi_dQ = dLambdaRi_dQ(:,ind);
                dLambda_dQ(j,:) = dLambdaRi_dQ(:).';
            otherwise
        end
    end

    % P
    if dims.nP
        % Block formulation
        switch method
            case 'rudisill'
                slv = C_rud\[zeros(1,dims.nP);(kron(vtf(thetaP,dims.nx,dims.nu),(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtf(thetaP,dims.nx,dims.nu),eye(r)))];
                dLambda_dP(j,:) = slv(end,:);
                dZj_dP2 = slv(1:end-1,:);
                dl_dP(ptr:ptr+m-1,:) = C*dZj_dP2 + Ls*kron(eye(m),(sys.Q*Xr(:,j)).');
                dXr_dP_rsh(:,j) = dZj_dP2(:);
            case 'murthy'
                slv = C_rud\(kron(vtf(thetaP,dims.nx,dims.nu),(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtf(thetaP,dims.nx,dims.nu),eye(r)));
                slv = [slv(1:m_max(j)-1,:);zeros(1,size(slv,2));slv(m_max(j):end,:)];
                dLambda_dP(j,:) = slv(end,:);
                dZj_dP2 = slv(1:end-1,:);
                dl_dP(ptr:ptr+m-1,:) = C*dZj_dP2 + Ls*kron(eye(m),(sys.Q*Xr(:,j)).');
                dXr_dP_rsh(:,j) = dZj_dP2(:);
            case 'magnus'
                dZj_dP = -pinv(lambda(j)*E-A,1e-12)*(eye(r,r)-(E*Xr(:,j)*Xl(:,j)')/(Xl(:,j)'*E*Xr(:,j)))*(kron(vtf(thetaP,dims.nx,dims.nu),(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtf(thetaP,dims.nx,dims.nu),eye(r)));
                dl_dP(ptr:ptr+m-1,:) = C*dZj_dP + Ls*kron(eye(m),(sys.Q*Xr(:,j)).');
                dXr_dP_rsh(:,j) = dZj_dP(:);

                dLambdaRi_dP = -Xl(:,j)'*(kron(vtf(thetaP,dims.nx,dims.nu),(sys.Q*Xr(:,j)).')+kron((sys.Q*Xr(:,j)).'*vtf(thetaP,dims.nx,dims.nu),eye(r)))/(Xl(:,j)'*E*Xr(:,j));
                dLambda_dP(j,:) = dLambdaRi_dP(:).';
            otherwise
        end
    end

    ptr = ptr + m;
end

% Compute (Gr-V*Gs.').'*Erinv.'*(dXrinv/dx).'*ej
if dims.nE, drR1_dE = -Xr\(reshape(permute(reshape((dXr_dE_rsh*EXinvB).',m,r,[]),[2,1,3]),r,[])); end
if dims.nJ, drR1_dJ = -Xr\(reshape(permute(reshape((dXr_dJ_rsh*EXinvB).',m,r,[]),[2,1,3]),r,[])); end
if dims.nR, drR1_dR = -Xr\(reshape(permute(reshape((dXr_dR_rsh*EXinvB).',m,r,[]),[2,1,3]),r,[])); end
if dims.nQ, drR1_dQ = -Xr\(reshape(permute(reshape((dXr_dQ_rsh*EXinvB).',m,r,[]),[2,1,3]),r,[])); end
if dims.nP, drR1_dP = -Xr\(reshape(permute(reshape((dXr_dP_rsh*EXinvB).',m,r,[]),[2,1,3]),r,[])); end

% Compute right residual derivativs dr_dx
pntr = 1;
for j=1:r
    if dims.nE
        slv = -B.'*Erinv.'*(kron(vtu(thetaE)',(EXinv(:,j)).')+kron((EXinv(:,j)).'*vtu(thetaE)',eye(r)));
        ind = nonzeros(tril(indMat));
        slv = slv(:,ind);
        dr_dE(pntr:pntr+m-1,:)=reshape(drR1_dE(j,:).',m,[])+slv;
    end
    if dims.nJ,dr_dJ(pntr:pntr+m-1,:)=reshape(drR1_dJ(j,:).',m,[]); end
    if dims.nR,dr_dR(pntr:pntr+m-1,:)=reshape(drR1_dR(j,:).',m,[]); end
    if dims.nQ,dr_dQ(pntr:pntr+m-1,:)=reshape(drR1_dQ(j,:).',m,[]); end
    if dims.nP,dr_dP(pntr:pntr+m-1,:)=reshape(drR1_dP(j,:).',m,[])-Ls*kron(eye(m),EXinv(:,j).'); end
    pntr=pntr+m;
end

% Bring to form dPhiR_dx = [dl1_dx;dr1_dx;dl2_dx;...;drr_dx]
ptr1 = 1;
ptr2 = 1;
for j=1:r
    if dims.nE,dPhi_dE(ptr1:ptr1+m-1,:) = dl_dE(ptr2:ptr2+m-1,:); dPhi_dE(ptr1+m:ptr1+2*m-1,:) = dr_dE(ptr2:ptr2+m-1,:); end
    if dims.nJ,dPhi_dJ(ptr1:ptr1+m-1,:) = dl_dJ(ptr2:ptr2+m-1,:); dPhi_dJ(ptr1+m:ptr1+2*m-1,:) = dr_dJ(ptr2:ptr2+m-1,:); end
    if dims.nR,dPhi_dR(ptr1:ptr1+m-1,:) = dl_dR(ptr2:ptr2+m-1,:); dPhi_dR(ptr1+m:ptr1+2*m-1,:) = dr_dR(ptr2:ptr2+m-1,:); end
    if dims.nQ,dPhi_dQ(ptr1:ptr1+m-1,:) = dl_dQ(ptr2:ptr2+m-1,:); dPhi_dQ(ptr1+m:ptr1+2*m-1,:) = dr_dQ(ptr2:ptr2+m-1,:); end
    if dims.nP,dPhi_dP(ptr1:ptr1+m-1,:) = dl_dP(ptr2:ptr2+m-1,:); dPhi_dP(ptr1+m:ptr1+2*m-1,:) = dr_dP(ptr2:ptr2+m-1,:); end

    ptr1 = ptr1+2*m;
    ptr2 = ptr2+m;
end

%% Derivatives w.r.t. matrices where dLambdaR_dx = 0
% G
if dims.nG
    ptr=1;
    for i=1:r
        dPhi_dG(ptr:(ptr+m-1),:)=kron(eye(m,m),(sys.Q*Xr(:,i)).');     % dlj_dG
        dPhi_dG((ptr+m):(ptr+2*m-1),:)=kron(eye(m,m),EXinv(:,i).'); % drj_dG
        ptr=ptr+2*m;
    end
end

% L
if dims.nL
    q = dims.nL/r;
    L = vtf(thetaL,r,q);
    I = eye(r);
    col=1;
    ptr=1;
    for k=1:q
        for j=1:r
            Ejk = zeros(r,q);Ejk(j,k)=1;
            dQr_dL = lyap(A',Ejk*L'+L*Ejk');
            for i=1:r
                dPhi_dL(ptr:(ptr+m-1),col) = sys.G'*dQr_dL*Xr*I(:,i);
                ptr = ptr + 2*m;
            end
            ptr=1;
            col=col+1;
        end
    end
end

%% Construct final matrices
dLambda_dx = [];
dPhi_dx = [];
if dims.nE, dLambda_dx = [dLambda_dx,dLambda_dE]; dPhi_dx = [dPhi_dx,dPhi_dE]; end
if dims.nJ, dLambda_dx = [dLambda_dx,dLambda_dJ]; dPhi_dx = [dPhi_dx,dPhi_dJ]; end
if dims.nR, dLambda_dx = [dLambda_dx,dLambda_dR]; dPhi_dx = [dPhi_dx,dPhi_dR]; end
if dims.nQ, dLambda_dx = [dLambda_dx,dLambda_dQ]; dPhi_dx = [dPhi_dx,dPhi_dQ]; end
if dims.nG, dLambda_dx = [dLambda_dx,dLambda_dG]; dPhi_dx = [dPhi_dx,dPhi_dG]; end
if dims.nP, dLambda_dx = [dLambda_dx,dLambda_dP]; dPhi_dx = [dPhi_dx,dPhi_dP]; end
if dims.nL, dLambda_dx = [dLambda_dx,dLambda_dL]; dPhi_dx = [dPhi_dx,dPhi_dL]; end

end

function [dLambda_dx,dPhi_dx,dl_dx,dr_dx,dXr_dx_rsh] = prealloc(r,m,dim)
dPhi_dx = zeros(2*r*m,dim);
dLambda_dx = zeros(r,dim);
dl_dx = zeros(r*m,dim);
dr_dx = zeros(r*m,dim);
dXr_dx_rsh = zeros(r*dim,r);
end
