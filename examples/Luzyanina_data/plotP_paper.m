% plotP plots the profiles stored in parameters.
%
% USAGE:
% ======
% fh = plotP(parameters)
% fh = plotP(parameters,fh)
% fh = plotP(parameters,fh,I)
% fh = plotP(parameters,fh,I,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%       and profiles.
% fh ... handle of figure in which profiles is plotted. If no
%       figure handle is provided, a new figure is opened.
% I ... index of subplot (parameter) which is updated. If no index is
%       provided the profiles for all parameters are updated.
% options ... options of plotting
%       .mark_constraint ... if 'true', points on the profile for which one
%           ore more parameters are close to the constraints are marked
%           with a cross (default = 'false').
%       .scale ... scale of y-axis (default = 'lin').
%
% Outputs:
% ========
% fh .. figure handle
%
% 2012/05/31 Jan Hasenauer

% function fh = plotP(parameters,fh,I,options)
function fh = plotP_paper(varargin)

%% CHECK AND ASSIGN INPUTS
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotPL requires a parameter object as input.');
end

% Open figure
if nargin >= 2
    if ~isempty(varargin{2})
        fh = figure(varargin{2});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Index of subplot which is updated
I = 1:length(parameters.P);
if nargin >= 3
    if ~isempty(varargin{3})
        I = varargin{3};
    end
end

% Options
options.mark_constraint = 'false';
options.alpha_level = 0.95;
options.alpha_level_type = 'single'; % 'simultanous'
options.type = 'R';
options.interval = 'dynamic';%'static';
options.hold_on = 'false';
options.linewidth = 1.5;
options.parameter_number = [];
options.plot_local = 'true';
if nargin == 4
    options = setdefault(varargin{4},options);
end

if ~isempty(options.parameter_number)
    npar_alpha = options.parameter_number;
else
    npar_alpha = parameters.number;
end

%% INITALIZATION
switch options.alpha_level_type
    case 'single'
        dof = 1;
    case 'simultanous'
        dof = parameters.number;
end

% Maximum a posterior estimate
logPost_max = max(parameters.MS.MAP_list.logPost);

%% CONSTRUCT INTERIOR BOUNDS
p = 1e-4;
xmin = (1-p)*parameters.min +    p *parameters.max;
xmax =    p *parameters.min + (1-p)*parameters.max;

%% PLOT PROFILE LIKELIHOODS
s = options.s;
lind = options.lind;

% Loop: Parameter
for l = 1:length(I)
    i = I(l);
    % Open subplot
    subplot(s(1),s(2),lind(l));
    if strcmp(options.hold_on,'true');
        hold on;
    else
        hold off;
    end
    
    % Plot profile likelihood
    if ~isempty(parameters.P(i).par)
        switch options.type
            case 'R'
                plot(parameters.P(i).par(i,:),exp(parameters.P(i).logPost - logPost_max),'-','linewidth',options.linewidth,'color',options.PL_color); hold on;
            case 'PL'
                plot(parameters.P(i).par(i,:),parameters.P(i).logPost,'-','linewidth',options.linewidth,'color',options.PL_color); hold on;
        end
    end
    
    % Determine index of points which are in the interior
    if strcmp(options.mark_constraint,'true')
        ind = find(sum(bsxfun(@gt,xmin,parameters.P(i).par)+bsxfun(@gt,parameters.P(i).par,xmax),1));
        if ~isempty(ind)
            switch options.type
                case 'R'
                    plot(parameters.P(i).par(i,ind),exp(parameters.P(i).logPost(ind) - logPost_max),'x','linewidth',options.linewidth,'color',options.threshold_color); hold on;    
                case 'PL'
                    plot(parameters.P(i).par(i,ind),parameters.P(i).logPost(ind),'x','linewidth',options.linewidth,'color',options.threshold_color); hold on;    
            end
        end
    end
    
    % Plot optimum
    ind = 1;%find(parameters.MS.MAP_list.logPost >= parameters.MS.MAP.logPost-chi2inv(options.alpha_level,dof)/2);
    switch options.type
        case 'R'
            plot(parameters.MS.MAP_list.par(i,ind),exp(parameters.MS.MAP_list.logPost(ind)-logPost_max),'o','linewidth',options.linewidth,'color',options.opt_color); hold on;
        case 'PL'
            plot(parameters.MS.MAP_list.par(i,ind),parameters.MS.MAP_list.logPost(ind),'o','linewidth',options.linewidth,'color',options.opt_color); hold on;
    end
    
    % Plot confidence levels

    % Limits
    switch options.interval
        case 'static'
            xl = [parameters.min(i),parameters.max(i)];
        case 'dynamic'
            if (size(parameters.P(i).par,2) >= 2) || (size(parameters.MS.MAP_list.par(i,ind),2) >= 2)
                xl = [min([parameters.P(i).par(i,:),min(parameters.MS.MAP_list.par(i,ind))]),...
                      max([parameters.P(i).par(i,:),max(parameters.MS.MAP_list.par(i,ind))])];
            else
                xl = [parameters.min(i),parameters.max(i)];
            end
        otherwise
            error('This option is not available.');
    end
    xlim(xl);
    switch options.type
        case 'R'
            ylim([0,1.1]);
        case 'PL'
            %ylim([parameters.MS.MAP.logPost-5,parameters.MS.MAP.logPost]);
    end            
    
    % Confidence levels
    switch options.type
        case 'R'
            plot(xl,[1,1]*exp(-chi2inv(options.alpha_level,1)/2),'k--');
            if strcmp(options.alpha_level_type,'simultanous')
                plot(xl,[1,1]*exp(-chi2inv(options.alpha_level,parameters.number)/2),'k:');
            end
        case 'PL'
            plot(xl,[1,1]*(parameters.MS.MAP.logPost-chi2inv(options.alpha_level,1)/2),'k--');
            if strcmp(options.alpha_level_type,'simultanous')
                plot(xl,[1,1]*(parameters.MS.MAP.logPost-chi2inv(options.alpha_level,parameters.number)/2),'k:');
            end
    end
    
    % Plot local approximation
    if strcmp(options.plot_local,'true')
        Sigma = inv(parameters.MS.MAP.hessian);
        sigma = sqrt(Sigma(i,i));
        x = parameters.P(i).par(i,:);
        x_opt = parameters.MS.MAP.par(i);
        switch options.type
            case 'R'
                plot(x,exp(-0.5*((x-x_opt)/sigma).^2),'b--','linewidth',options.linewidth); hold on;
            case 'PL'
                plot(x,parameters.MS.MAP.logPost-0.5*((x-x_opt)/sigma).^2,'b--','linewidth',options.linewidth); hold on;
        end
    end
    
    % Labels
    xlabel(parameters.name(i),'fontsize',10);
%     if mod(lind(l),s(2)) == 1 || lind(l) == 31
        switch options.type
            case 'R'
                ylabel({'likelihood','ratio, R'});
            case 'PL'
                ylabel('log-profile, log(PL)');
        end
%     else
%         set(gca,'Ytick',[]);
%     end
    set(gca,'FontSize',6);

end

drawnow;
