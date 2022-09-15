function [p,v] = samPassive(sys,varargin)
% samPassive - Check if LTI system sys is passive or not by sampling along
%              the imaginary axis and checking positive definiteness of the
%              Popov function.
%
% Syntax:
%   p = samPassive(sys)
%   p = samPassive(sys,w)
%   p = samPassive(sys,opts)
%   p = samPassive(sys,w,opts)
%
% Description:
%       p = samPassive(sys,w) evaluates the Popov function at the frequencies
%       in w along the imaginary axis and checks if all eigenvalues are greater
%       than zero.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      sss object, containing LTI system
%       *Optional Input Arguments:*
%       - w:    Vector containing sampling points. 
%               [{logspace(-15,5,5000)} / vector]
%       - Opts: Structure with execution parameters
%           - .tol:             Tolerance for passivity verification. If,
%                               for example, Psi has eigenvalues close to
%                               zero, it will still be considered positive
%                               if its absolute value is less than Opts.tol.
%                               [{1e-7} / double]
%           - .plot:            Defines if a plot of the Popov function is created
%                               after passivity verification.
%                               [{false} / true]
%           - .outputVolume:    Defines the volume of the output regarding
%                               passivity violations. 'bands': Passivity
%                               violating bands, 'all': every sampling
%                               point, 'subset': every k-th sampling
%                               point,'worst': only the maximum violation.
%                               [{'bands'} / 'all' / 'subset' / 'worst'].
%           - .subFreq:         Defines the output volume if 'subset' is
%                               selected for Opts.outputVolume.
%                               [{4} / integer]
%
% Output Arguments:
%       - p:    Flag specifiying if system is passive or not
%       - v:    Structure containing passivity violation informations
%
% References:
%       [1] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: 
%           Theory and Applications, 2016.
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

%% Input parsing
[sys,w,opts] = parseInputs(sys,varargin{:});

%% Evaluate Popov function
% Preallocation
lambda = zeros(1,length(w));
nue = zeros(size(sys.B,2),length(w));

% Evaluation
for i = 1:length(w)
    % Transfer function
    G=evalG(sys,1i*w(i));

    % Popov function
    Psi=G+G';

    % Collect eigenvalues and corresponding eigenvectors
    [V,D] = eig(full(Psi));
    [lambda(i),idx] = min(real(diag(D)));
    nue(:,i) = V(:,idx);
end

%% Check passivity
% Lift Popov function by tolerance
x = lambda+opts.tol;
p=true;

% Check if violations are present
if any(x<0)
    p=false;
end

% Form output
if strcmp(opts.outputVolume, 'bands')
    % Get all sign changes
    idx_pn = find(lambda(1:end-1)>=0 & lambda(2:end) < 0);
    idx_np = find(lambda(1:end-1)<0 & lambda(2:end) >= 0);

    % Add up all sign changes and end points
    idx = [idx_pn,idx_np];
    idx = sort(idx,'ascend');
    idx_1 = find(idx==1,1);
    if isempty(idx_1)
        idx = [1,idx,numel(w)];
    end

    % Preallocation
    bands = zeros(numel(idx)-1,2);
    vl = zeros(numel(idx)-1,1);
    f = zeros(numel(idx)-1,1);
    ev = zeros(size(sys.B,2),numel(idx)-1);

    % Define bands, worst case violations and corresponding frequencies
    for i=1:numel(idx)-1
        bands(i,1) = w(idx(i));
        bands(i,2) = w(idx(i+1));
        vl(i) = min(lambda(idx(i)+1:idx(i+1)));
        k = find(lambda==vl(i),1);
        f(i) = w(k);
        ev(:,i) = nue(:,k);
    end

    % Define output
    v.bands = bands;
    v.violations = vl;
    v.frequencies = f;
    v.eigenvectors = ev;

elseif strcmp(opts.outputVolume, 'all')
    % Return all information
    v.violations = lambda;
    v.frequencies = w;
    v.eigenvectors = nue;

elseif strcmp(opts.outputVolume, 'subset')
    omega = w;

    % Get worst violations
    [w_lambda,idx] = mink(lambda,30);
    w_omega = omega(idx);
    w_nue = nue(:,idx);

    lambda(idx) = [];
    omega(idx) = [];
    nue(:,idx) = [];

    z = 1:opts.subFreq:length(omega);
    v.violations = [lambda(z),w_lambda];
    v.frequencies = [omega(z),w_omega];
    v.eigenvectors = [nue(:,z),w_nue];

elseif strcmp(opts.outputVolume, 'worst')
    % Get worst violations
    [wc_lambda,idx] = min(lambda);
    wc_w = w(idx);
    wc_nue = nue(:,idx);

    v.violations = wc_lambda;
    v.frequencies = wc_w;
    v.eigenvectors = wc_nue;

end

%% Plot lambda curve
if opts.plot
    fg=figure();
    fg.Position = [200 200 900 600];

    hold on
    grid on

    % Filter negative eigenvalues
    xbelow = x<0;

    % Define colors
    green = [0 1 .2];
    red = [1 0 .2];

    % Plot point in (0,0) for legend entries
    plot(0,0,'-','Color',green,'LineWidth', 2, 'MarkerSize', 2)
    plot(0,0,'-','Color',red,'LineWidth', 2, 'MarkerSize', 2)

    % Plot data
    plot(w(~xbelow),x(~xbelow),'.','Color',green,'LineWidth', 2, 'MarkerSize', 2)
    plot(w(xbelow),x(xbelow),'.','Color',red,'LineWidth', 2, 'MarkerSize', 2)

    % Set figure and axis properties
    %     title('Popov function $\Psi(j \omega)$','interpreter','latex')
    xlabel('$\omega$','Interpreter','latex')
    ylabel('$\lambda_{min}(\Psi(j \omega))$','Interpreter','latex')
    set(gca, 'XScale', 'log')
    set(gca,'TickLabelInterpreter','latex')
    legend('$\lambda_{min}(\Psi(j \omega))>0$','$\lambda_{min}(\Psi(j \omega))<0$','Interpreter','Latex','Location','northoutside')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    xticks(logspace(-20,20,41));
    xlim([w(1) w(end)])
    ax=gca;
    ax.YMinorGrid = 'off';
    ax.XMinorGrid = 'off';
end

end

%% Auxilliary functions
function G_si = evalG(sys,s_i)
% Evaluates the transfer function G(s) at certain shifts s_i.
res = (s_i*sys.E-sys.A)\sys.B;
G_si = sys.C*res+sys.D;
end

function [sys,w,opts] = parseInputs(sys,varargin)
if nargin < 2
    opts = struct();
    w = logspace(-15,5,5000);
    w = [0,w];
elseif nargin < 3 && ~isa(varargin{1},'struct')
    opts = struct();
    w = varargin{1};
elseif nargin < 3
    opts = varargin{1};
    w = logspace(-15,5,5000);
    w = [0,w];
else
    if isa(varargin{1},'struct')
        opts = varargin{1};
        w = varargin{2};
    else
        opts = varargin{2};
        w = varargin{1};
    end
end

% Admissible options set
OptsAdmissible.tol = 1e-7;
OptsAdmissible.plot = {true,false};
OptsAdmissible.outputVolume = {'bands','all','subset','worst'};
OptsAdmissible.subFreq = 4;
opts = phsMOR_parseOpts(opts, OptsAdmissible);

if isa(sys,'phs') || isa(sys,'phsRed')
    warning('System is port-Hamiltonian, hence passive.')
    sys = phs2sss(sys);
end

% Make sure E exists
if isempty(sys.E)
    sys.E = eye(size(sys.A));
end

end
