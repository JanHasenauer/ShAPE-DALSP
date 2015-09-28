% getMultiStarts.m computes the maximum a posterior estimate of the
%   parameters of a user-supplied posterior function. Therefore, a
%   multi-start local optimization is used.
%
% Note: This function can exploit up to (n_start + 1) workers when running
% in 'parallel' mode.
%
% USAGE:
% ======
% [...] = getMultiStarts(parameters,objective_function)
% [...] = getMultiStarts(parameters,objective_function,options)
% [parameters,fh] = getMultiStarts(...)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least:
%   .number ... number of parameter
%   .guess ... initial guess of parameter
%   .min ... lower bound for parameter values
%   .max ... upper bound for parameter values
%   .name = {'name1',...} ... names of the parameters
%   .init_fun ... function to draw starting points for local
%   	optimization. The function has to have the input structure
%           .init_fun(theta_0,theta_min,theta_max)
%       Alternatively, a latin hypercube or a uniform random sampling can
%       be used by setting the respective options
% objective_function ... objective function to be optimized. This function
%       should possess exactly one input, the parameter vector.
% options ... options of algorithm
%   .obj_type ... type of objective function provided
%       = 'log-posterior' (default) ... algorithm assumes that
%               log-posterior or log-likelihood are provided and perfroms
%               a maximization of the objective function.
%       = 'negative log-posterior' ... algorithm assumes that negative
%               log-posterior or negative log-likelihood are provided and
%               perfroms a minimization of the objective function.
%   .comp_type ... type of computations
%       = 'sequential' (default) ... classical sequential (in core) method
%       = 'parallel' ... multi-core method exploiting parfor
%   .fmincon ... options for fmincon (the local optimizer)
%   .n_starts ... number of local optimizations (default = 20).
%   .init_threshold ... log-likelihood / log-posterior threshold for 
%       initialization of optimization (default = -inf).
%   .proposal ... method used to propose starting points
%       = 'latin hypercube' (default) ... latin hypercube sampling
%       = 'uniform' ... uniform random sampling
%       = 'user-supplied' ... user supplied function parameters.init_fun
%   .rng ... initialization of random number generator (default = 0).
%       = any ral number r => random generator is initialized with r.
%       = [] ... random number generator is not initialized.
%       (Initializing the random number generator can be helpfult to
%       understand problems.)
%   .mode ... output of algorithm
%       = 'visual' (default) ... plots are gnerated which show the progress
%       = 'text' ... optimization results for multi-start is printed on screen
%       = 'silent' ... no output during the multi-start local optimization
%   .fh ... handle of figure in which results are printed. If no
%       handle is provided, a new figure is used.
%   .plot_options ... plot options for plotMultiStarts.m.
%   .save ... determine whether results are directly saved
%       = false (default) ... results are not saved
%       = true ... results are stored do an extra folder
%   .trace ... determine whether objective function and parameter values
%       are stored over iterations
%       = false (default) ...  not saved
%       = true ... stored in .par_trace and .fval_trace
%   .foldername ... name of the folder in which results are stored.
%       If no folder is provided, a random foldername is generated.
%
% Outputs:
% ========
% parameters ... updated parameter object containing:
%   .MS ... information about multi-start optimization
%       .par(:,i) ... ith MAP
%       .par0(:,i) ... starting point yielding ith MAP
%       .logPost(i) ... log-posterior for ith MAP
%       .logPost0(i) ... log-posterior for starting point yielding ith MAP
%       .gradient(:,i) ... gradient of log-posterior at ith MAP
%       .hessian(:,:,i) ... hessian of log-posterior at ith MAP
%       .n_objfun(i) ... # objective evaluations used to calculate ith MAP
%       .n_iter(i) ... # iterations used to calculate ith MAP
%       .t_cpu(i) ... CPU time for calculation of ith MAP
%       .exitflag(i) ... exitflag the optimizer returned for ith MAP
%       .par_trace(:,:,i) ... parameter trace for ith MAP
%       .fval_trace(:,i) ... objective function value trace for ith MAP
% fh ... figure handle
%
% 2012/05/31 Jan Hasenauer
% 2012/07/11 Jan Hasenauer
% 2014/06/11 Jan Hasenauer
% 2015/07/28 Fabian Froehlich

% function [parameters,fh] = getMultiStarts(parameters,objective_function,options)
function [parameters,fh] = getMultiStarts(varargin)

%% Check inputs and assign default values
if nargin >= 2
    parameters = varargin{1};
    objective_function = varargin{2};
else
    error('getMultiStarts.m requires at least two inputs.')
end

% Check parameters:
if ~isfield(parameters,'min') || ~isfield(parameters,'max')
    error('Algorithm requires lower and upper bounds');
else
    parameters.min = parameters.min(:);
    parameters.max = parameters.max(:);
end
if length(parameters.min) ~= length(parameters.max)
    error('Dimension of parameters.min and parameters.max does not agree.');
else
    if max(parameters.min >= parameters.max)
        error('There exists at least one i for which parameters.min(i) >= parameters.max(i).');
    end
end
if ~isfield(parameters,'number')
    parameters.number = length(parameters.min);
else
    if parameters.number ~= length(parameters.min)
        error('Dimension mismatch: parameters.number ~= length(parameters.min).');
    end
end
if isfield(parameters,'guess')
    if ~isempty(parameters.guess);
        if size(parameters.guess,1) ~= length(parameters.max)
            error('Dimension of parameters.guess does not agree with dimesion of parameters.min and .max.');
        end
    end
end
constr.A = [];
constr.b = [];
constr.Aeq = [];
constr.beq = [];
if isfield(parameters,'constraints')
    parameters.constraints = setdefault(parameters.constraints,constr);
else
    parameters.constraints = constr;
end

% Check initial guess
if ~isfield(parameters,'guess')
    parameters.guess = [];
end

% Check and assign options
options.fmincon = optimset('algorithm','interior-point',...
    'display','off',...
    'GradObj','on',...
    'PrecondBandWidth',inf);
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.n_starts = 20; % number of multi-start local optimizations
options.proposal = 'latin hypercube'; % 'uniform','user-supplied'
options.init_threshold = -inf;
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.rng = 0;
options.foldername = strrep(datestr(now,31),' ','__');
options.fh = [];
options.trace = false;
options.plot_options.add_points.par = [];
if nargin == 3
    options = setdefault(varargin{3},options);
end


%% Initialization and figure generation
fh = [];
switch options.mode
    case 'visual'
        if isempty(options.fh)
            fh = figure;
        else
            fh = figure(options.fh);
        end
    case 'text'
        fprintf(' \nOptimization:\n=============\n');
    case 'silent' % no output
        % Force fmincon to be silent.
        options.fmincon = optimset(options.fmincon,'display','off');
end

%% Initialization of random number generator
if ~isempty(options.rng)
    rng(options.rng);
end

%% Sampling of starting points
switch options.proposal
    case 'latin hypercube'
        % Sampling from latin hypercube
        par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                             lhsdesign(options.n_starts - size(parameters.guess,2),parameters.number,'smooth','off')'))];
    case 'uniform'
        % Sampling from latin hypercube
        par0 = [parameters.guess,...
                bsxfun(@plus,parameters.min,bsxfun(@times,parameters.max - parameters.min,...
                             rand(parameters.number,options.n_starts - size(parameters.guess,2))))];
    case 'user-supplied'
        % Sampling from user-supplied function
        par0 = [parameters.guess,...
                parameters.init_fun(parameters.guess,parameters.min,parameters.max,options.n_starts - size(parameters.guess,2))];
end
parameters.MS.n_starts = options.n_starts;
parameters.MS.par0 = par0;

%% Initialization
parameters.MS.par = nan(parameters.number,options.n_starts);
parameters.MS.logPost0 = nan(options.n_starts,1);
parameters.MS.logPost = nan(options.n_starts,1);
parameters.MS.gradient = nan(parameters.number,options.n_starts);
parameters.MS.hessian  = nan(parameters.number,parameters.number,options.n_starts);
parameters.MS.n_objfun = nan(options.n_starts,1);
parameters.MS.n_iter = nan(options.n_starts,1);
parameters.MS.t_cpu = nan(options.n_starts,1);
parameters.MS.exitflag = nan(options.n_starts,1);
if(options.trace)
    parameters.MS.par_trace = nan(parameters.number,options.fmincon.MaxIter+1,options.n_starts);
    parameters.MS.fval_trace = nan(options.fmincon.MaxIter+1,options.n_starts);
end

%% Preperation of folder
if options.save
    try
        rmdir(options.foldername,'s');
    end
    mkdir(options.foldername);
    save([options.foldername '/init'],'parameters');
end

%% Multi-start local optimization -- SEQUENTIAL
if strcmp(options.comp_type,'sequential')
    
    global error_count
    % initialise tracing of parameter and objective function values
    if(options.trace)
        options.fmincon.OutputFcn = @outfun_fmincon_wtrace;
    else
        options.fmincon.OutputFcn = @outfun_fmincon_wotrace;
    end
    
    % Loop: Mutli-starts
    for i = 1:options.n_starts
        % Reset error count
        error_count = 0;
        
        % Evaluation of objective function at starting point
        if strcmp(options.fmincon.GradObj,'off')
            J_0 = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
        elseif strcmp(options.fmincon.GradObj,'on') && ~strcmp(options.fmincon.Hessian,'user-supplied')
            [J_0,grad_J_0] = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
        else
            [J_0,grad_J_0,H_J_0] = obj_w_error_count(parameters.MS.par0(:,i),objective_function,options.obj_type);
        end
        parameters.MS.logPost0(i) = -J_0;
        
        % Optimization
        t_cpu_fmincon = cputime;
        if J_0 < -options.init_threshold
            if(options.trace)
                history.x = NaN(parameters.number,options.fmincon.MaxIter+1);
                history.fval = NaN(options.fmincon.MaxIter+1,1);
            end
            
            % Optimization using fmincon
            [theta,J_opt,parameters.MS.exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
                fmincon(@(theta) obj_w_error_count(theta,objective_function,options.obj_type),...  % negative log-likelihood function
                par0(:,i),...    % initial parameter
                parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                parameters.min,...     % lower bound
                parameters.max,...     % upper bound
                [],options.fmincon);   % options
            
            % Assignment
            parameters.MS.logPost(i) = -J_opt;
            parameters.MS.par(:,i) = theta;
            parameters.MS.gradient(:,i) = gradient_opt;
            if isempty(hessian_opt)
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                end
            elseif max(hessian_opt(:)) == 0
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                end
            end
            parameters.MS.hessian(:,:,i) = full(hessian_opt);
            parameters.MS.n_objfun(i) = results_fmincon.funcCount;
            parameters.MS.n_iter(i) = results_fmincon.iterations;
            if(options.trace)
                parameters.MS.par_trace(:,:,i) = history.x;
                parameters.MS.fval_trace(:,i) = history.fval;
            end
        end
        parameters.MS.t_cpu(i) = cputime - t_cpu_fmincon;
        
        % Save
        if options.save
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__logPost.csv'],parameters.MS.logPost(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__logPost0.csv'],parameters.MS.logPost0(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__par.csv'],parameters.MS.par(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__par0.csv'],parameters.MS.par0(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__gradient.csv'],parameters.MS.gradient(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__hessian.csv'],parameters.MS.hessian(:,:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__t_cpu.csv'],parameters.MS.t_cpu(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__n_objfun.csv'],parameters.MS.n_objfun(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__n_iter.csv'],parameters.MS.n_iter(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__exitflag.csv'],parameters.MS.exitflag(i),'delimiter',',','precision',12);
        end
        
        % Output
        switch options.mode
            case 'visual', fh = plotMultiStarts(parameters,fh,options.plot_options);
            case 'text', disp(['  ' num2str(i,'%d') '/' num2str(options.n_starts,'%d')]);
            case 'silent' % no output
        end
    end
    
    % Assignment
    parameters = sortMultiStarts(parameters);
    
end

%% Multi-start local optimization -- PARALLEL
if strcmp(options.comp_type,'parallel')
    
    % Initialization
    par = nan(parameters.number,options.n_starts);
    logPost0 = nan(options.n_starts,1);
    logPost = nan(options.n_starts,1);
    gradient = nan(parameters.number,options.n_starts);
    hessian  = nan(parameters.number,parameters.number,options.n_starts);
    n_objfun = nan(options.n_starts,1);
    n_iter = nan(options.n_starts,1);
    t_cpu = nan(options.n_starts,1);
    exitflag = nan(options.n_starts,1);
    
    % Loop: Mutli-starts
    parfor i = 1:options.n_starts
        % Evaluation of objective function at starting point
        if strcmp(options.fmincon.GradObj,'off')
            J_0 = obj(par0(:,i),objective_function,options.obj_type);
        elseif strcmp(options.fmincon.GradObj,'on') && ~strcmp(options.fmincon.Hessian,'user-supplied')
            [J_0,grad_J_0] = obj(par0(:,i),objective_function,options.obj_type);
        else
            [J_0,grad_J_0,H_J_0] = obj(par0(:,i),objective_function,options.obj_type);
        end
        logPost0(i) = -J_0;
        
        % Optimization
        t_cpu_fmincon = cputime;
        if J_0 < -options.init_threshold
            % Optimization using fmincon
            [theta,J_opt,exitflag(i),results_fmincon,~,gradient_opt,hessian_opt] = ...
                fmincon(@(theta) obj(theta,objective_function,options.obj_type),...  % negative log-posterior function
                par0(:,i),...    % initial parameter
                parameters.constraints.A  ,parameters.constraints.b  ,... % linear inequality constraints
                parameters.constraints.Aeq,parameters.constraints.beq,... % linear equality constraints
                parameters.min,...     % lower bound
                parameters.max,...     % upper bound
                [],options.fmincon);   % options
            
            % Assignment
            logPost(i) = -J_opt;
            par(:,i) = theta;
            gradient(:,i) = gradient_opt;
            if isempty(hessian_opt)
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                end
            elseif max(abs(hessian_opt(:))) == 0
                if strcmp(options.fmincon.Hessian,'user-supplied')
                    [~,~,hessian_opt] = obj(theta,objective_function,options.obj_type);
                end
            end
            hessian(:,:,i) = full(hessian_opt);
            n_objfun(i) = results_fmincon.funcCount;
            n_iter(i) = results_fmincon.iterations;
        end
        t_cpu(i) = cputime - t_cpu_fmincon;
        
        % Save
        if options.save
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__logPost.csv'],logPost(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__logPost0.csv'],logPost0(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__par.csv'],par(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__par0.csv'],par0(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__gradient.csv'],gradient(:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__hessian.csv'],hessian(:,:,i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__t_cpu.csv'],t_cpu(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__n_objfun.csv'],n_objfun(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__n_iter.csv'],n_iter(i),'delimiter',',','precision',12);
            dlmwrite([options.foldername '/MS' num2str(i,'%d') '__exitflag.csv'],exitflag(i),'delimiter',',','precision',12);
        end
        
        % Output
        switch options.mode
            case 'text', disp(['  ' num2str(i,'%d') '/' num2str(options.n_starts,'%d')]);
            case {'silent','visual'} % no output
        end
    end
    
    % Assignment
    parameters.MS.par0 = par0;
    parameters.MS.par = par;
    parameters.MS.logPost0 = logPost0;
    parameters.MS.logPost = logPost;
    parameters.MS.gradient = gradient;
    parameters.MS.hessian  = hessian;
    parameters.MS.n_objfun = n_objfun;
    parameters.MS.n_iter = n_iter;
    parameters.MS.t_cpu = t_cpu;
    parameters.MS.exitflag = exitflag;
    parameters = sortMultiStarts(parameters);
    
    % Output
    switch options.mode
        case 'visual', fh = plotMultiStarts(parameters,fh);
        case {'text','silent'} % no output
    end
    
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Multi-start optimization FINISHED.');
    case 'silent' % no output
end

%% Nested function for storing of objective function and parameter values
    function stop = outfun_fmincon_wtrace(x,optimValues,state)
        
        switch state
            case 'init'
                % do nothing
            case 'interrupt'
                % do nothing
            case 'iter'
                if optimValues.iteration>0
                    history.fval(optimValues.iteration) = optimValues.fval;
                    history.x(:,optimValues.iteration) = x;
                end
            case 'done'
                % do nothing
        end
        
        if error_count <= 20
            stop = false;
        else
            warning('Too many failed objective function evaluations.')
            stop = true;
        end
        
    end

    function stop = outfun_fmincon_wotrace(x,optimValues,state)
        
        switch state
            case 'init'
                % do nothing
            case 'interrupt'
                % do nothing
            case 'iter'
                % do nothing
            case 'done'
                % do nothing
        end
        
        if error_count <= 20
            stop = false;
        else
            warning('Too many failed objective function evaluations.')
            stop = true;
        end
        
    end
end


%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = obj(theta,fun,type)

try
    switch nargout
        case {0,1}
            J = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J,-G};
                case 'negative log-posterior' , varargout = { J, G};
            end
        case 3
            [J,G,H] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J,-G,-H};
                case 'negative log-posterior' , varargout = { J, G, H};
            end
            if(any(isnan(H)))
                error('Hessian contains NaNs')
            end
    end
    
catch error_msg
    % disp(error_msg.message)
    
    % Derive output
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(theta),1)};
        case 3
            varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = obj_w_error_count(theta,fun,type)

global error_count

try
    switch nargout
        case {0,1}
            J = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J,-G};
                case 'negative log-posterior' , varargout = { J, G};
            end
        case 3
            [J,G,H] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {-J,-G,-H};
                case 'negative log-posterior' , varargout = { J, G, H};
            end
            if(any(isnan(H)))
                error('Hessian contains NaNs')
            end
    end
    % Reset error count
    error_count = 0;
catch error_msg
    % Increase error count
    error_count = error_count + 1;
    
    % Display a warning with error message
    warning(['Evaluation of likelihood failed because: ' error_msg.message]);
    
    % Derive output
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(theta),1)};
        case 3
            varargout = {inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

