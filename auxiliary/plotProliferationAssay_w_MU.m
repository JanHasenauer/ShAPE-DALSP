% plotProliferationAssay plots the outcome of a proliferation assay
%     experiment.
%
% USAGE:
% ======
% fh = plotProliferationAssay_w_MU(S,D)
% fh = plotProliferationAssay_w_MU(S,D,fh)
% fh = plotProliferationAssay_w_MU(S,D,fh,options)
%
% INPUTS:
% =======
% S ... simulation result of DSLP or DALSP model, containg at least .Mhist.
%       If S is empty, only the experimental data are plotted.
% D ... data object.
%       If D is empty, only the simulation is plotted.
% fh ... handle of figure in which profile likelihood is plotted. If no
%       figure handle is provided, a new figure is opened.
% options ... options of plotting
%       .lincol ... colors of line.
%
% Outputs:
% ========
% fh .. figure handle
%
% 2012/06/04 Jan Hasenauer

% function fh = plotProliferationAssay_w_MU(S,D,fh,options)
function fh = plotProliferationAssay_w_MU(varargin)

%% CHECK AND ASSIGN INPUTS
% Assign parameters
if nargin >= 2
    S = varargin{1};
    D = varargin{2};
else
    error('plotProliferationAssay requires a simuation and a data object as input.');
end
% Check if inputs are empty
if isempty(S)*isempty(D)
    error('S or D must be non-empty')
else
    kmax = 0;
    if ~isempty(S)
        kmax = size(S.n_y,1);
    end
    if ~isempty(D)
        kmax = max(kmax,length(D.t_plot));
    end
end

% Open figure
if nargin >= 3
    if ~isempty(varargin{3})
        fh = figure(varargin{3});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Options
%options.lincol = {[0 0 0],[0 0 1],[1 0 0],[0,0.8,0],[0 1 1]};
options.lincol = colormap(jet(kmax));
options.plot_simulation = 'true';
options.plot_data = 'true';
options.label1.x = 'CFSE fluorescence [UI]';
options.label1.y = 'cell counts';
options.label2.x = 'CFSE fluorescence [UI]';
options.label2.y = 'p([CFSE flu.]|time)*[CFSE flu.]';
options.label3.x = 'time';
options.label3.y = 'number of cells';
options.holdoff = 'true';
options.two_plots = 'true';
if nargin == 4
    options = setdefault(varargin{4},options);
end

if isempty(S)
    options.plot_simulation = 'false';
end
if isempty(D)
    options.plot_data = 'false';
end

%% INITIALIZATION
lh = [];
ls = {};
if strcmp(options.plot_data,'true')
    x = [D(1).bins(1,1); D(1).bins(:,2)]';
else
    x = S.x;
end
if strcmp(options.plot_simulation,'true')
    if ~isfield(S,'H')
        S.H = bsxfun(@times,S.n_y(:,1:end-1) + S.n_y(:,2:end), 0.5*diff(x(:)'));
    end
    S.h = bsxfun(@times,S.p_y(:,1:end-1) + S.p_y(:,2:end), 0.5*diff(x(:)'));
    lsd = 'o';
    lwd = 1;
else
    lwd = 1;
    lsd = '-o';
end

% Calculation of lower and upper bounds
lb = zeros(size(S.H));
ub = zeros(size(S.H));
for j = 1:length(D.t)
    for i = 1:length(lb)
        b = binoinv([0.05 0.95],sum(D.cellcount(j,:)),S.h(j,i));
        lb(j,i) = b(1);
        ub(j,i) = b(2);
    end
end

%% CFSE DISTRIBUTION
if strcmp(options.two_plots,'false')
if strcmp(options.holdoff,'true')
    hold off;
end
% Plot mueasurement and simulation / construct legend label
for j = 1:length(D.t)
    % Measurement
    if strcmp(options.plot_data,'true')
        lh(end+1) = stairs(x,[D.cellcount(j,:),0],...
                  '-','Color',min(0.2+options.lincol(mod(j-1,length(options.lincol))+1,:),1),...
                  'linewidth',lwd); hold on;
        ls{end+1} = [D.t_name{j} ' - measurement'];
    end
    % Simulation
    if strcmp(options.plot_simulation,'true')
        lh(end+1)=stairs(x,S.H(j,[1:end,end]),...
                  '-','Color',0.8*options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',2);
        stairs(x,lb(j,[1:end,end]),...
                  '--','Color',0.8*options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',1);
        stairs(x,ub(j,[1:end,end]),...
                  '--','Color',0.8*options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',1);
        ls{end+1} = [D.t_name{j} ' - simulation'];
    end
end

% Set legend / label
set(gca,'Xscale','log');
xlim(x([1,end]));
xlabel(options.label1.x);
ylabel(options.label1.y);
legend(lh,ls,'location','Best');
end

%% CFSE DISTRIBUTION
if strcmp(options.two_plots,'true')
subplot(1,3,1:2);
if strcmp(options.holdoff,'true')
    hold off;
end
% Plot mueasurement and simulation / construct legend label
for j = 1:length(D.t)
%     % Measurement
%     if strcmp(options.plot_data,'true')
%         lh(end+1) = stairs(x,(D.cellcount(j,[1:end,end])/sum(D.cellcount(j,:))),...
%                   '-','Color',0.7*options.lincol(mod(j-1,length(options.lincol))+1,:),...
%                   'linewidth',lwd); hold on;
%         ls{end+1} = [D.t_name{j} ' - measurement'];
%     end
%     % Simulation
%     if strcmp(options.plot_simulation,'true')
%         lh(end+1)=stairs(x,S.h(j,[1:end,end]),...
%                   '-','Color',options.lincol(mod(j-1,length(options.lincol))+1,:),...
%                   'linewidth',2);
% 
%         indx = round(1:0.5:length(x)-0.5);
%         indb = round(0.5:0.5:length(x)-1);
%         fill(x([indx,indx(end:-1:1)]),[lb(j,indb),ub(j,indb(end:-1:1))]/sum(D.cellcount(j,:)),...
%                   min(0.5+options.lincol(mod(j-1,length(options.lincol))+1,:),1),'facealpha',.2,...
%                   'edgecolor',options.lincol(mod(j-1,length(options.lincol))+1,:),'edgealpha',.2);
%         %set(h,'facealpha',.5)

    % Measurement
    if strcmp(options.plot_data,'true')
        lh(end+1) = stairs(x,(D.cellcount(j,[1:end,end])/sum(D.cellcount(j,:))),...
                  '-','Color',min(0.2+options.lincol(mod(j-1,length(options.lincol))+1,:),1),...
                  'linewidth',lwd); hold on;
        ls{end+1} = [D.t_name{j} ' - measurement'];
    end
    % Simulation
    if strcmp(options.plot_simulation,'true')
        lh(end+1)=stairs(x,S.h(j,[1:end,end]),...
                  '-','Color',options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',2);
        stairs(x,lb(j,[1:end,end])/sum(D.cellcount(j,:)),...
                  '--','Color',0.8*options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',1);
        stairs(x,ub(j,[1:end,end])/sum(D.cellcount(j,:)),...
                  '--','Color',0.8*options.lincol(mod(j-1,length(options.lincol))+1,:),...
                  'linewidth',1);
        ls{end+1} = [D.t_name{j} ' - simulation'];
    end
end

% Set legend / label
set(gca,'Xscale','log');
xlim(x([1,end]));
xlabel(options.label2.x);
ylabel(options.label2.y);
legend(lh,ls,'location','Best');

%% CELL NUMBER
subplot(1,3,3);
if strcmp(options.holdoff,'true')
    hold off;
end
% Plot mueasurement and simulation
ls = {};
if strcmp(options.plot_simulation,'true')
    plot(D.t+D.t_plot(1),S.N,'-','linewidth',2,...
        'Color',0.8*options.lincol(mod(0,length(options.lincol))+1,:)); hold on;
    ls{end+1} = 'simulation';
end
if strcmp(options.plot_data,'true')
    plot(D.t_plot,sum(D.cellcount,2),lsd,'linewidth',2,...
        'Color',min(0.2+options.lincol(mod(0,length(options.lincol))+1,:),1)); hold on;
    ls{end+1} = 'measurement';
end

% Set legend / label
xlim(D.t_plot([1,end]));
xlabel(options.label3.x);
ylabel(options.label3.y);
set(gca,'yscale','log');
if length(ls) == 2
    legend(ls,'location','Best');
end
end


drawnow;
