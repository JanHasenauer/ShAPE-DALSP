% plotParameterUncertainty.m visualizes profile likelihood and MCMC samples
% stored in parameters.
%
% USAGE:
% ======
% fh = plotParameterUncertainty(parameters,type)
% fh = plotParameterUncertainty(parameters,type,fh)
% fh = plotParameterUncertainty(parameters,type,fh,I)
% fh = plotParameterUncertainty(parameters,type,fh,I,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%       and results of optimization (.MS) and uncertainty analysis
%       (.P and .S). This structures is the output of plotMultiStarts.m,
%       getProfiles.m or plotSamples.m.
% type ... string indicating the type of visualization: '1D' or '2D'
% fh ... handle of figure. If no figure handle is provided, a new figure
%       is opened.
% I ... index of parameters which are updated. If no index is provided
%       all parameters are updated.
% options ... options of plotting
%   .hold_on ... indicates whether plots are redrawn or whether something
%       is added to the plot
%       = 'false' (default) ... new plot
%       = 'true' ... extension of plot
%   .interval ... selection mechanism for x limits
%       = 'dynamic' (default) ... x limits depending on analysis results
%       = 'static' ... x limits depending on parameters.min and .max or on
%          user-defined bound options.bounds.min and .max. The later are
%          used if provided.
%   .bounds ... bounds used for visualization if options.interval = 'static'
%       .min ... lower bound
%       .max ... upper bound
%   .P ... options for profile plots
%       .plot_type ... plot type
%           = 0 (default if no profiles are provided) ... no plot
%           = 1 (default if profiles are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of profile lines (default = [1,0,0])
%       .lw ... line width of profile lines (default = 2)
%       .name ... name of legend entry (default = 'P')
%   .S ... options for sample plots
%       .plot_type ... plot type
%           = 0 (default if no samples are provided) ... no plot
%           = 1 (default if samples are provided) ... histogram
%           = 2 ... kernel-density estimates
%       .hist_col ... color of histogram (default = [0.7,0.7,0.7])
%       .bins ... number of histogram bins
%           = 'optimal' ... selection using Scott's rule
%           = 'conservative' (default) ... selection using Scott's rule / 2
%           = N (with N being an integer) ... N bins
%       .sp_col ... color of scatter plot (default = [0.7,0.7,0.7])
%       .sp_m ... marker for scatter plot (default = '.')
%       .sp_ms ... marker size for scatter plot (default = 5)
%       .name ... name of legend entry (default = 'S')
%   .MS ... options for multi-start optimization plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio and
%               position of optima above threshold
%           = 2 ... negative log-likelihood and position of optima 
%               above threshold
%       .col ... color of local optima (default = [1,0,0])
%       .lw ... line width of local optima (default = 1.5)
%       .name_conv ... name of legend entry (default = 'MS - conv.')
%       .name_nconv ... name of legend entry (default = 'MS - not conv.')
%       .only_optimum ... only optimum is plotted
%   .A ... options for distribution approximation plots
%       .plot_type ... plot type
%           = 0 (default if no MS are provided) ... no plot
%           = 1 (default if MS are provided) ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .col ... color of approximation lines (default = [0,0,1])
%       .lw ... line width of approximation lines (default = 2)
%       .sigma_level ... sigma-level which is visualized (default = 2)
%       .name ... name of legend entry (default = 'P_{app}')
%   .boundary ... options for boundary visualization
%       .mark ... marking of profile points which are on the boundary
%           = 0 ... no visualization
%           = 1 (default) ... indicates points which ar close to the
%               boundaries in one or more dimensions.
%       .eps ... minimal distance from boundary for which points are
%               consider to e close do the boundary (default = 1e-4). Note
%               that a one-norm is used.
%   .CL ... options for confidence level plots
%       .plot_type ... plot type
%           = 0 (default) ... no plot
%           = 1 ... likelihood ratio
%           = 2 ... negative log-likelihood
%       .alpha ... visualized confidence level (default = 0.95)
%       .type ... type of confidence interval
%           = 'point-wise' (default) ... point-wise confidence interval
%           = 'simultanous' ... point-wise confidence interval
%           = {'point-wise','simultanous'} ... both
%       .col ... color of profile lines (default = [0,0,0])
%       .lw ... line width of profile lines (default = 2)
%       .name ... name of legend entry (default = 'cut-off'):
%   .op2D ... options used for 2D plot to position subplot axes.
%       .b1 ... offset from left and bottom border (default = 0.15)
%       .b2 ... offset from left and bottom border (default = 0.02)
%       .r ... relative width of subplots (default = 0.95)
%   .add_points ... option used to add additional points, e.g. true
%           parameter in the case of test examples
%       .par ... n x m matrix of m additional points
%       .col ... color used for additional points (default = [0,0.8,0]).
%                  This can also be a m x 3 matrix of colors.
%       .ls ... line style (default = '--')
%       .lw ... line width (default = 2)
%       .m ... marker style (default = 'd')
%       .ms ... line width (default = 8)
%       .name ... name of legend entry (default = 'add. point')
%   .legend ... legend options
%       .color ... background color (default = 'none').
%       .box ... legend outine (default = 'on').
%       .orientation ... orientation of list (default = 'vertical').
%   .fontsize ... fontsize
%       .tick ... fontsize for ticklabels (default = 12).
%
% Outputs:
% ========
% fh ... figure handle
%
% 2012/05/31 Jan Hasenauer
% 2014/06/20 Jan Hasenauer

% function fh = plotParameterUncertainty(parameters,type,fh,I,options)
function fh = plotParameterUncertainty(varargin)

%% Check and assign inputs
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotParameterUncertainty requires a parameter object as input.');
end

% Plot type
type = '1D';
if nargin >= 2
    if ~isempty(varargin{2})
        type = varargin{2};
    end
end
if ~max(strcmp({'1D','2D'},type))
    error('The ''type'' of plot is unknown.')
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

% Index of subplot which is updated
I = 1:parameters.number;
if nargin >= 4
    if ~isempty(varargin{4})
        I = varargin{4};
        if ~isnumeric(I) || max(abs(I - round(I)) > 0)
            error('I is not an integer vector.');
        end
    end
end

% Options
% General plot options
options.hold_on = 'false';
options.interval = 'dynamic'; %'static';

% Default profile plotting options
%   0 => no plot
%   1 => likelihood ratio
%   2 => negative log-likelihood
if isfield(parameters,'P')
    options.P.plot_type = 1;
else
    options.P.plot_type = 0; 
end
options.P.col = [1,0,0];
options.P.lw = 2;
options.P.name = 'P';

% Default sample plotting options
%   0 => no plot
%   1 => histogram
%   2 => kernel-density estimate
if isfield(parameters,'S')
    options.S.plot_type = 1; 
else
    options.S.plot_type = 0; 
end
options.S.bins = 'conservative';
options.S.scaling = [];
options.S.hist_col = [0.7,0.7,0.7];
options.S.sp_col = [0.7,0.7,0.7];
options.S.lin_col = [1,0,0];
options.S.lin_lw = 2;
options.S.sp_m = '.';
options.S.sp_ms = 5;
options.S.PT.lw = 1.5;
options.S.PT.ind = [];
options.S.PT.col = [];
options.S.PT.plot_type = 0;
if isfield(parameters,'S')
    if isfield(parameters.S,'PT');
        options.S.PT.plot_type = options.S.plot_type;
        options.S.PT.ind = 1:size(parameters.S.PT.par,3);
        options.S.PT.col = [linspace(0,1,size(parameters.S.PT.par,3))',...
                            0.2*ones(size(parameters.S.PT.par,3),1),...
                            linspace(1,0,size(parameters.S.PT.par,3))'];
    end
end
options.S.name = 'S';

% Local optima
%   0 => no plot
%   1 => likelihood ratio
%   2 => negative log-likelihood
if isfield(parameters,'MS')
    if options.S.plot_type == options.P.plot_type
        options.MS.plot_type = options.S.plot_type;
    elseif options.S.plot_type == 0
        options.MS.plot_type = options.P.plot_type;
    elseif options.P.plot_type == 0
        options.MS.plot_type = options.S.plot_type;
    end
else
    options.MS.plot_type = 0; 
end
options.MS.col = [1,0,0];
options.MS.lw = 2;
options.MS.name_conv = 'MS - conv.';
options.MS.name_nconv = 'MS - not conv.';
options.MS.only_optimum = false;

% Default approxiamtion plotting options
%   0 => no plot
%   1 => likelihood ratio
%   2 => negative log-likelihood
options.A.plot_type = options.MS.plot_type;
options.A.col = [0,0,1];
options.A.lw = 2;
options.A.sigma_level = 2;
options.A.name = 'P_{app}';

% Boundary detection
options.boundary.mark = 1;
options.boundary.eps = 1e-4;

% Confidence level
options.CL.plot_type = 0;%options.MS.plot_type;
options.CL.alpha = 0.95;
options.CL.type = 'point-wise'; % 'simultanous', {'point-wise','simultanous'}
options.CL.col = [0,0,0];
options.CL.lw = 2;
options.CL.name = 'cut-off';

% Settings for 2D plot
options.op2D.b1 = 0.15;
options.op2D.b2 = 0.02;
options.op2D.r = 0.95;

% Additional points
options.add_points.par = [];
options.add_points.col = [0,0.8,0];
options.add_points.ls = '--';
options.add_points.lw = 2;
options.add_points.m = 's';
options.add_points.ms = 8;
options.add_points.name = 'add. point';

% Legend
options.legend.color = 'none';
options.legend.box = 'on';
options.legend.orientation = 'vertical';
options.legend.position = [];

% Labels
options.labels.y_always = true;
options.labels.y_name = [];

% Fontsize
options.fontsize.tick = 12;

% Assignment of user-provided options
if nargin == 5
    options = setdefault(varargin{5},options);
end

% Check
if ~isfield(parameters,'P')
    options.boundary.mark = 0;
end

% Subplot arrangement
if ~isfield(options,'subplot_size_1D')
    options.subplot_size_1D = round(sqrt(length(I))*[1,1]);
    if prod(options.subplot_size_1D) < length(I)
        options.subplot_size_1D(2) = options.subplot_size_1D(2) + 1;
    end
end
if ~isfield(options,'subplot_indexing_1D')
    options.subplot_indexing_1D = 1:length(I);
end

%% INITALIZATION
% Maximum a posterior estimate
if isfield(parameters,'MS')
    logPost_max = max(parameters.MS.logPost);
end

% Degrees of freedom (for chi^2 test)
dof = 1;
if max(strcmp(options.CL.type,'simultanous'))
    dof = parameters.number;
end

%% 1D Parameter distributions
if strcmp(type,'1D')

% Compute number of subfigure

% Loop: Parameter
for l = 1:length(I)
    % Initialization of legend
    legh = [];
    legs = {};

    % Assign parameter index
    i = I(l);
    
    % Open subplot
    subplot(options.subplot_size_1D(1),options.subplot_size_1D(2),options.subplot_indexing_1D(l));
    
    % Hold on/off
    if strcmp(options.hold_on,'true')
        hold on;
    else
        hold off;
    end
    
    % Boundaries
    switch options.interval
        case 'dynamic'
            xl = [+inf,-inf];
            
            if isfield(parameters,'MS')
                if max(strcmp(options.CL.type,'point-wise'))
                    L = find(parameters.MS.logPost(:) > (parameters.MS.logPost(1)-chi2inv(options.CL.alpha,1)/2));
                end
                if max(strcmp(options.CL.type,'simultanous'))
                    L = find(parameters.MS.logPost(:) > (parameters.MS.logPost(1)-chi2inv(options.CL.alpha,parameters.numbe)/2));
                end
                xl(1) = min(xl(1),min(parameters.MS.par(i,L)));
                xl(2) = max(xl(2),max(parameters.MS.par(i,L)));
            end
        
            flag_plot_P = 0;
            if options.P.plot_type >= 1
                if length(parameters.P) >= i
                    if ~isempty(parameters.P(i).par)
                        xl(1) = min(xl(1),min(parameters.P(i).par(i,:)));
                        xl(2) = max(xl(2),max(parameters.P(i).par(i,:)));
                        flag_plot_P = 1;
                    end
                end
            end

            if options.S.plot_type >= 1
                xl(1) = min(xl(1),min(parameters.S.par(i,:)));
                xl(2) = max(xl(2),max(parameters.S.par(i,:)));
            end
            
            if xl(1) == xl(2)
                xl(1) = xl(1) - 1e-10;
                xl(2) = xl(2) + 1e-10;
            end
        case 'static'
            if isfield(options,bounds)
                xl = [options.bounds.min(i),options.bounds.max(i)];
            else
                xl = [parameters.min(i),parameters.max(i)];
            end
    end

    % Plot: Visualizaion of MCMC samples of tempered posterior distribution
    switch options.S.PT.plot_type
        case 0
            % no plot
        case 1
            % histogram
            if isfield(parameters.S,'PT') && options.S.PT.plot_type
                for k = options.S.PT.ind
                    switch options.S.bins
                        case 'optimal'
                            h = 3.49*std(parameters.S.PT.par(i,:,k))/(length(parameters.S.PT.par(i,:,k))^(1/3));
                            nbin = round((max(parameters.S.PT.par(i,:,k))-min(parameters.S.PT.par(i,:,k)))/h);
                        case 'conservative'
                            h = 2*3.49*std(parameters.S.PT.par(i,:,k))/(length(parameters.S.PT.par(i,:,k))^(1/3));
                            nbin = round((max(parameters.S.PT.par(i,:,k))-min(parameters.S.PT.par(i,:,k)))/h);
                        otherwise
                            nbin = options.S.bins;
                    end
                    [N,X] = hist(parameters.S.PT.par(i,:,k),nbin);
                    bar(X,N/max(N),1,'facecolor','none','edgecolor',options.S.PT.col(k,:)); hold on;
                end
            end
        case 2
            % kernel-density estimate
            if isfield(parameters.S,'PT') && options.S.PT.plot_type
                for k = options.S.PT.ind
                    x_grid = linspace(min(parameters.S.PT.par(i,:,k)),max(parameters.S.PT.par(i,:,k)),100);
                    [KDest] = kde_simple(squeeze(parameters.S.PT.par(i,:,k)),x_grid);
                    plot(x_grid,KDest/max(KDest),'-','color',options.S.PT.col(k,:),'linewidth',options.S.PT.lw); hold on;
                end
            end
        otherwise
            error('Selected value for ''options.S.plot_type'' is not available.');
    end
        
    % Plot: Visualizaion of MCMC samples of posterior distribution
    h = [];
    switch options.S.plot_type
        case 0
            % no plot
        case 1
            % histogram
            switch options.S.bins
                case 'optimal'
                    b = 3.49*std(parameters.S.par(i,:))/(length(parameters.S.par(i,:))^(1/3));
                    nbin = round((max(parameters.S.par(i,:))-min(parameters.S.par(i,:)))/b);
                case 'conservative'
                    b = 2*3.49*std(parameters.S.par(i,:))/(length(parameters.S.par(i,:))^(1/3));
                    nbin = round((max(parameters.S.par(i,:))-min(parameters.S.par(i,:)))/b);
                otherwise
                    nbin = options.S.bins;
            end
            [N,X] = hist(parameters.S.par(i,:),nbin);
            h = bar(X,N/max(N),1,'facecolor',options.S.hist_col); hold on;
        case 2
            % kernel-density estimate
            x_grid = linspace(min(parameters.S.par(i,:)),max(parameters.S.par(i,:)),100);
            [KDest] = kde_simple(squeeze(parameters.S.par(i,:)),x_grid);
            h = plot(x_grid,KDest/max(KDest),'-','color',options.S.lin_col,'linewidth',options.S.lin_lw); hold on;
        otherwise
            error('Selected value for ''options.S.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.S.name;
    end

    % Plot: Local approximation
    h = [];
    switch options.A.plot_type
        case 0
            % no plot
        case 1
            if isfield(parameters.MS,'hessian')
                % likelihood ratio
                Sigma = pinv(parameters.MS.hessian(:,:,1));
                sigma = sqrt(Sigma(i,i));
                % Get grid
                par_grid = parameters.MS.par(i,1) + sigma*linspace(-4,4,100);
                par_grid = par_grid(find((parameters.min(i) <= par_grid).*(par_grid <= parameters.max(i))));
                % Calculation of objectiev function approximation
                % - with non-zero gradient
%                 ind_I = [1:i-1,i+1:parameters.number];
%                 dtheta_i = -parameters.MS.par(i,1)+par_grid;
%                 dtheta_ind_I = -pinv(parameters.MS.hessian(ind_I,ind_I,1))*bsxfun(@plus,parameters.MS.hessian(ind_I,i,1)*dtheta_i,parameters.MS.gradient(ind_I,1));
%                 dtheta = [dtheta_ind_I(1:i-1,:);dtheta_i;dtheta_ind_I(i:end,:)];
%                 for l = 1:size(dtheta,2)
%                     dtheta(:,l) = max(min(parameters.MS.par(:,1)+dtheta(:,l),parameters.max),parameters.min)-parameters.MS.par(:,1);
%                 end
%                 J = nan(1,size(dtheta,2));
%                 for l = 1:size(dtheta,2)
%                     J(l) = parameters.MS.gradient(:,1)'*dtheta(:,l) + 0.5*dtheta(:,l)'*parameters.MS.hessian(:,:,1)*dtheta(:,l);
%                 end
                J = parameters.MS.gradient(i,1)*(par_grid-parameters.MS.par(i,1)) + 0.5*((par_grid-parameters.MS.par(i,1))/sigma).^2;
                % - with zero gradient
%                 J = 0.5*((par_grid-parameters.MS.par(i,1))/sigma).^2;
                % Plot
                h = plot(par_grid,exp(-J),'-','linewidth',options.A.lw,'color',options.A.col); hold on;
            else
                warning('No hessian provided in .MS. Approximation in not plotted.');
            end
        case 2
            if isfield(parameters.MS,'hessian')
                % negative log-likelihood
                Sigma = pinv(parameters.MS.hessian(:,:,1));
                sigma = sqrt(Sigma(i,i));
                % Get grid
                par_grid = parameters.MS.par(i,1) + sigma*linspace(-4,4,100);
                par_grid = par_grid(find((parameters.min(i) <= par_grid).*(par_grid <= parameters.max(i))));
                % Calculation of objectiev function approximation
                % - with non-zero gradient
%                 ind_I = [1:i-1,i+1:parameters.number];
%                 dtheta_i = -parameters.MS.par(i,1)+par_grid;
%                 dtheta_ind_I = -pinv(parameters.MS.hessian(ind_I,ind_I,1))*bsxfun(@plus,parameters.MS.hessian(ind_I,i,1)*dtheta_i,parameters.MS.gradient(ind_I,1));
%                 dtheta = [dtheta_ind_I(1:i-1,:);dtheta_i;dtheta_ind_I(i:end,:)];
%                 J = nan(1,size(dtheta,2));
%                 for l = 1:size(dtheta,2)
%                     J(l) = parameters.MS.gradient(:,1)'*dtheta(:,l) + 0.5*dtheta(:,l)'*parameters.MS.hessian(:,:,1)*dtheta(:,l);
%                 end
                J = parameters.MS.gradient(i,1)*(par_grid-parameters.MS.par(i,1)) + 0.5*((par_grid-parameters.MS.par(i,1))/sigma).^2;
                % - with zero gradient
%                 J = 0.5*((par_grid-parameters.MS.par(i,1))/sigma).^2;
                % Plot
                h = plot(par_grid,J,'-','linewidth',options.A.lw,'color',options.A.col); hold on;
            else
                warning('No hessian provided in .MS. Approximation in not plotted.');
            end
        otherwise
            error('Selected value for ''options.A.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.A.name;
    end

    % Plot: Profile likelihood
    h = [];
    switch options.P.plot_type * flag_plot_P
        case 0
            % no plot
        case 1
            % likelihood ratio
            h = plot(parameters.P(i).par(i,:),exp(parameters.P(i).logPost - logPost_max),'-','linewidth',options.P.lw,'color',options.P.col); hold on;
        case 2
            % negative log-likelihood
            h = plot(parameters.P(i).par(i,:),parameters.P(i).logPost,'-','linewidth',options.P.lw,'color',options.P.col); hold on;
        otherwise
            error('Selected value for ''options.P.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.P.name;
    end
    
    % Plot: Additional points
    h = [];
    if ~isempty(options.add_points.par)
        % Check dimension:
        if size(options.add_points.par,1) ~= parameters.number
            warning(['The matrix options.add_points.par should possess ' num2str(parameters.number) ' rows.']);
        else
            for j = 1:size(options.add_points.par,2)
                if size(options.add_points.col,1) == size(options.add_points.par,2)
                    l = j;
                else
                    l = 1;
                end
                h = plot(options.add_points.par(i,j)*[1,1],[0,1.05],options.add_points.ls,'color',options.add_points.col(l,:),'linewidth',options.add_points.lw);
            end
        end
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.add_points.name;
    end
    
    % Bounds
    if (options.P.plot_type >= 1) * flag_plot_P
    switch options.boundary.mark
        case 0
            % no plot
        case 1
            ind = find(sum( bsxfun(@gt,parameters.min+options.boundary.eps,parameters.P(i).par)...
                           +bsxfun(@gt,parameters.P(i).par,parameters.max-options.boundary.eps),1));
            if ~isempty(ind)
                switch options.P.plot_type
                    case 1
                        % likelihood ratio
                        plot(parameters.P(i).par(i,ind),exp(parameters.P(i).logPost(ind) - logPost_max),'x','linewidth',options.P.lw,'color',options.P.col); hold on;    
                    case 2
                        % negative log-likelihood
                        plot(parameters.P(i).par(i,ind),parameters.P(i).logPost(ind),'x','linewidth',options.P.lw,'color',options.P.col); hold on;    
                end
            end
        otherwise
            error('Selected value for ''options.boundary.mark'' is not available.');
    end
    end
    
    % Plot: Local optima
    h_conv = [];
    h_nconv = [];
    if options.MS.only_optimum
        ind = 1;
    else
        ind = find(parameters.MS.logPost >= parameters.MS.logPost(1)-chi2inv(options.CL.alpha,dof)/2);
    end
    ind_conv = ind(find(min((parameters.MS.exitflag(ind) > 0)+(parameters.MS.exitflag(ind) == -3),1)));
    ind_nconv = setdiff(ind,ind_conv);
    switch options.MS.plot_type
        case 0
            % no plot
        case 1
            % likelihood ratio
            h_conv = plot(parameters.MS.par(i,ind_conv),exp(parameters.MS.logPost(ind_conv)-logPost_max),'o','linewidth',options.MS.lw,'color',options.MS.col); hold on;
            h_nconv = plot(parameters.MS.par(i,ind_nconv),exp(parameters.MS.logPost(ind_nconv)-logPost_max),'s','linewidth',options.MS.lw,'color',options.MS.col);
        case 2
            % negative log-likelihood
            h_conv = plot(parameters.MS.par(i,ind_conv),parameters.MS.logPost(ind_conv),'o','linewidth',options.MS.lw,'color',options.MS.col); hold on;
            h_nconv = plot(parameters.MS.par(i,ind_nconv),parameters.MS.logPost(ind_nconv),'s','linewidth',options.MS.lw,'color',options.MS.col); hold on;
        otherwise
            error('Selected value for ''options.MS.plot_type'' is not available.');
    end
    if ~isempty(h_conv)
        legh(end+1) = h_conv;
        legs{end+1} = options.MS.name_conv;
    end
    if ~isempty(h_nconv)
        legh(end+1) = h_nconv;
        legs{end+1} = options.MS.name_nconv;
    end
    
    % Limits
    % x
    if strcmp(options.interval,'static')
        xl = [parameters.min(i),parameters.max(i)];
    end
    xlim(xl);

    % y
    switch options.P.plot_type
        case {0,1}
            % likelihood ratio
            ylim([0,1.1]);
        case 2
            % Best choice not clear => automatic assignment
    end
    
    % Plot: Confidence levels
    h = [];
    switch options.CL.plot_type
        case 0
            % no plot
        case 1
            % likelihood ratio
            if max(strcmp(options.CL.type,'point-wise'))
                plot(xl,[1,1]*exp(-chi2inv(options.CL.alpha,1)/2),'--','color',options.CL.col);
            end
            if max(strcmp(options.CL.type,'simultanous'))
                plot(xl,[1,1]*exp(-chi2inv(options.CL.alpha,parameters.number)/2),':','linewidth',options.CL.lw,'color',options.CL.col);
            end
        case 2
            % negative log-likelihood
            if max(strcmp(options.CL.type,'point-wise'))
                plot(xl,[1,1]*(parameters.MS.logPost(1)-chi2inv(options.CL.alpha,1)/2),'--','linewidth',options.CL.lw,'color',options.CL.col);
            end
            if max(strcmp(options.CL.type,'simultanous'))
                plot(xl,[1,1]*(parameters.MS.logPost(1)-chi2inv(options.CL.alpha,parameters.number)/2),':','linewidth',options.CL.lw,'color',options.CL.col);
            end
        otherwise
            error('Selected value for ''options.CL.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.CL.name;
    end
    
    % Labels
    xlabel(parameters.name(i));
    if (mod(options.subplot_indexing_1D(l),options.subplot_size_1D(2)) == 1) || (length(I) == 1) || options.labels.y_always
        if isempty(options.labels.y_name)
            switch options.CL.plot_type
                case 0
                    % no plot
                    ylabel('post. prob., p');
                case 1
                    % likelihood ratio
                    ylabel('ratio, R');
                case 2
                    % negative log-likelihood
                    ylabel('log-profile, log(PL)');
            end
        else
            ylabel(options.labels.y_name);
        end
    else
        set(gca,'Ytick',[]);
    end
    set(gca,'fontsize',options.fontsize.tick);

    % Legend
    if l == 1
        if isempty(options.legend.position)
            legend(legh,legs,'color',options.legend.color,'box',options.legend.box,'orientation',options.legend.orientation);
        else
            legend(legh,legs,'color',options.legend.color,'box',options.legend.box,'orientation',options.legend.orientation,'position',options.legend.position);
        end
    end
end

end


%% 2D Parameter distributions
if strcmp(type,'2D')
    
% Loop: Parameter
for l1 = 1:length(I)
for l2 = 1:length(I)
    % Initialization of legend
    legh = [];
    legs = {};

    % Assign parameter index
    i1 = I(l1);
    i2 = I(l2);
    
    % Open subplot
%    subplot(length(I),length(I),(i2-1)*length(I)+i1);
    d = (1-options.op2D.b1-options.op2D.b2)/length(I);
    subplot('Position',[options.op2D.b1+(l1-1)*d,...
                        options.op2D.b1+(length(I)-l2)*d,...
                        options.op2D.r*d,options.op2D.r*d]);
    
    % Hold on/off
    if strcmp(options.hold_on,'true')
        hold on;
    else
        hold off;
    end
    
    % Boundaries
    switch options.interval
        case 'dynamic'
            xl1 = [+inf,-inf];
            xl2 = [+inf,-inf];
            
            flag_plot_P_i1 = 0;
            flag_plot_P_i2 = 0;
            if options.P.plot_type >= 1
                if length(parameters.P) >= i1
                    if ~isempty(parameters.P(i1).par)
                        xl1(1) = min(xl1(1),min(parameters.P(i1).par(i1,:)));
                        xl1(2) = max(xl1(2),max(parameters.P(i1).par(i1,:)));
                        flag_plot_P_i1 = 1;
                    end
                end
                if length(parameters.P) >= i2
                    if ~isempty(parameters.P(i2).par)
                        xl2(1) = min(xl2(1),min(parameters.P(i2).par(i2,:)));
                        xl2(2) = max(xl2(2),max(parameters.P(i2).par(i2,:)));
                        flag_plot_P_i2 = 1;
                    end
                end
            end

            if options.S.plot_type >= 1
                xl1(1) = min(xl1(1),min(parameters.S.par(i1,:)));
                xl1(2) = max(xl1(2),max(parameters.S.par(i1,:)));
                xl2(1) = min(xl2(1),min(parameters.S.par(i2,:)));
                xl2(2) = max(xl2(2),max(parameters.S.par(i2,:)));
            end

        case 'static'
            if isfield(options,bounds)
                xl1 = [options.bounds.min(i1),options.bounds.max(i1)];
                xl2 = [options.bounds.min(i2),options.bounds.max(i2)];
            else
                xl1 = [parameters.min(i1),parameters.max(i1)];
                xl2 = [parameters.min(i2),parameters.max(i2)];
            end
    end

    % Plot: MCMC samples
    h = [];
    switch options.S.plot_type
        case 0
            % no plot
        case 1
            % scatter plot
            h = plot(parameters.S.par(i1,:),parameters.S.par(i2,:),options.S.sp_m,...
                'color',options.S.sp_col,'markersize',options.S.sp_ms); hold on;
        otherwise
            error('Selected value for ''options.S.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.S.name;
    end

    % Plot: Local approximation
    h = [];
    switch options.A.plot_type
        case 0
            % no plot
        case {1,2}
            if isfield(parameters.MS,'hessian')
                Sigma = pinv(parameters.MS.hessian(:,:,1));
                % plot
                X = getEllipse(parameters.MS.par([i1,i2],1),Sigma([i1,i2],[i1,i2]),options.A.sigma_level);
                h = plot(X(1,:),X(2,:),'-','linewidth',options.A.lw/1.5,'color',options.A.col); hold on;
            else
                warning('No hessian provided in .MS. Approximation in not plotted.');
            end
        otherwise
            error('Selected value for ''options.A.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.A.name;
    end
    
    % Plot: Local optima
    h_conv = [];
    h_nconv = [];
    if options.MS.only_optimum
        ind = 1;
    else
        ind = find(parameters.MS.logPost >= parameters.MS.logPost(1)-chi2inv(options.CL.alpha,dof)/2);
    end
    ind_conv = ind(find(min((parameters.MS.exitflag(ind) > 0)+(parameters.MS.exitflag(ind) == -3),1)));
    ind_nconv = setdiff(ind,ind_conv);
    switch options.P.plot_type
        case 0
            % no plot
        case {1,2}
            h_conv = plot(parameters.MS.par(i1,ind_conv),parameters.MS.par(i2,ind_conv),'o','linewidth',options.MS.lw,'color',options.MS.col); hold on;
            h_nconv = plot(parameters.MS.par(i1,ind_nconv),parameters.MS.par(i2,ind_nconv),'s','linewidth',options.MS.lw,'color',options.MS.col); hold on;
        otherwise
            error('Selected value for ''options.MS.plot_type'' is not available.');
    end
    if ~isempty(h_conv)
        legh(end+1) = h_conv;
        legs{end+1} = options.MS.name_conv;
    end
    if ~isempty(h_nconv)
        legh(end+1) = h_nconv;
        legs{end+1} = options.MS.name_nconv;
    end

    % Plot: Profile likelihood
    h = [];
    switch options.P.plot_type
        case 0
            % no plot
        case {1,2}
            if flag_plot_P_i1
                h = plot(parameters.P(i1).par(i1,:),parameters.P(i1).par(i2,:),'-','linewidth',options.P.lw,'color',options.P.col*0.8); hold on;
            end
            if flag_plot_P_i2
                h = plot(parameters.P(i2).par(i1,:),parameters.P(i2).par(i2,:),'-','linewidth',options.P.lw,'color',options.P.col*0.6); hold on;
            end
        otherwise
            error('Selected value for ''options.P.plot_type'' is not available.');
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.P.name;
    end
    
    % Plot: Additional points
    h = [];
    if ~isempty(options.add_points.par)
        % Check dimension:
        if size(options.add_points.par,1) ~= parameters.number
            warning(['The matrix options.add_points.par should possess ' num2str(parameters.number) ' rows.']);
        else
            for j = 1:size(options.add_points.par,2)
                if size(options.add_points.col,1) == size(options.add_points.par,2)
                    l = j;
                else
                    l = 1;
                end
                h = plot(options.add_points.par(i1,j),options.add_points.par(i2,j),options.add_points.m,...
                    'color',options.add_points.col(l,:),'linewidth',options.add_points.lw,'markersize',options.add_points.ms);
            end
        end
    end
    if ~isempty(h)
        legh(end+1) = h;
        legs{end+1} = options.add_points.name;
    end
    
    % Bounds
    switch options.boundary.mark
        case 0
            % no plot
        case 1
            % i1
            if length(parameters.P) >= i1
                if ~isempty(parameters.P(i1).par)
                    ind = find(sum( bsxfun(@gt,parameters.min+options.boundary.eps,parameters.P(i1).par)...
                                    +bsxfun(@gt,parameters.P(i1).par,parameters.max-options.boundary.eps),1));
                    if ~isempty(ind)
                        switch options.P.plot_type
                            case {1,2}
                                plot(parameters.P(i1).par(i1,ind),parameters.P(i1).par(i2,ind),'x','linewidth',options.P.lw,'color',options.P.col*0.8); hold on;    
                        end
                    end
                end
            end
            
            % i2
            if length(parameters.P) >= i2
                if ~isempty(parameters.P(i2).par)
                    ind = find(sum( bsxfun(@gt,parameters.min+options.boundary.eps,parameters.P(i2).par)...
                                    +bsxfun(@gt,parameters.P(i2).par,parameters.max-options.boundary.eps),1));
                    if ~isempty(ind)
                        switch options.P.plot_type
                            case {1,2}
                                plot(parameters.P(i2).par(i1,ind),parameters.P(i2).par(i2,ind),'x','linewidth',options.P.lw,'color',options.P.col*0.6); hold on;    
                        end
                    end
                end
            end
        otherwise
            error('Selected value for ''options.boundary.mark'' is not available.');
    end
        
    % Limits
    if ~isinf(xl1(1))
        xlim(xl1);
    end
    if ~isinf(xl2(1))
        ylim(xl2);
    end

    % Labels
    if l2 == length(I)
        xlabel(parameters.name(i1));
    else
        set(gca,'xticklabel',[]);
    end
    if i1 == 1
        ylabel(parameters.name(i2));
    else
        set(gca,'yticklabel',[]);
    end
    set(gca,'fontsize',options.fontsize.tick);
    
    % Legend
    if (l1 == 1) && (l2 == 1)
        if isempty(options.legend.position)
            legend(legh,legs,'color',options.legend.color,'box',options.legend.box,'orientation',options.legend.orientation);
        else
            legend(legh,legs,'color',options.legend.color,'box',options.legend.box,'orientation',options.legend.orientation,'position',options.legend.position);
        end
    end

end
end

end


%% Update plot
drawnow;


