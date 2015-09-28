% function boxplot_horizontal(S,name,options)
function boxplot_horizontal(varargin)

%% Check inputs and assign default values
if nargin >= 2
    S = varargin{1};
    name = varargin{2};
else
    error('boxplot_horizontal.m requires at least two inputs.')
end

% Check parameters
if length(S) ~= length(name)
    error('''S'' and ''name'' have the same length.');
end

% Check and assign options
options.fill = true;
options.box_width = 0.3;
options.whisker_width = 0.15;
options.line_color = [0,0,0];
options.fill_color = [0.8,0.8,0.8];
options.line_width = 1;
options.line_width_median = 1;
options.marker_size = 6;
options.marker = 'o';
options.xscale = 'lin';
if nargin == 3
    options = setdefault(varargin{3},options);
end

% Loop: Samples
for i = 1:length(S)
    % Fill color
    if size(options.fill_color,1) == 1
        fc = options.fill_color;
    else
        fc = options.fill_color(i,:);
    end
    
    % Line color
    if size(options.line_color,1) == 1
        lc = options.line_color;
    else
        lc = options.line_color(i,:);
    end
    
    % Quantiles
    switch options.xscale
        case 'lin'
            Q = prctile(S{i},[25,50,75]);
            IQR = Q(3)-Q(1);
        case 'log'
            Q = prctile(log10(S{i}),[25,50,75]);
            IQR = Q(3)-Q(1);
            Whisker_min = min(S{i}(log10(S{i}) > (Q(1)-1.5*IQR)));
            Whisker_max = max(S{i}(log10(S{i}) < (Q(3)+1.5*IQR)));
            Q = 10.^Q;
    end
    
    % Whisker
    plot([Whisker_min,Whisker_max],length(S)-i*[1,1],'-','color',lc,'linewidth',options.line_width); hold on;
    plot(Whisker_min*[1,1],length(S)-i+options.whisker_width*[-1,1],'-','color',lc,'linewidth',options.line_width);
    plot(Whisker_max*[1,1],length(S)-i+options.whisker_width*[-1,1],'-','color',lc,'linewidth',options.line_width);
    
    % Box
    if options.fill
        fill([Q(1),Q(3),Q(3),Q(1)],length(S)-i+options.box_width*[-1,-1,1,1],lc,'facecolor',fc,'linewidth',options.line_width,'edgecolor',lc);
    else
        fill([Q(1),Q(3),Q(3),Q(1)],length(S)-i+options.box_width*[-1,-1,1,1],lc,'facecolor',[1,1,1],'linewidth',options.line_width,'edgecolor',lc);
    end

    % Median
    plot(Q(2)*[1,1],length(S)-i+options.box_width*[-1,1],'-','color',lc,'linewidth',options.line_width_median);
    
    % Outliers
    outliers = [S{i}(S{i} < Whisker_min);S{i}(Whisker_max < S{i})];
    if ~isempty(outliers)
        plot(outliers,length(S)-i,options.marker,'color',lc,'linewidth',options.line_width,'markersize',options.marker_size);
    end
end

% Label
set(gca,'ytick',[0:(length(S)-1)]);
if isempty(name)
    set(gca,'yticklabel',[length(S):-1:1]);
else
    set(gca,'yticklabel',name(end:-1:1));
end