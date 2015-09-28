% getParameterProfiles.m calculates the profiles of a user-supplied function,
%   starting from the maximum a posteriori estimate.
%
% Note: This function can exploit up to (n_theta + 1) workers when running
% in 'parallel' mode.
%
% USAGE:
% ======
% [...] = getParameterProfiles(parameters,objective_function)
% [...] = getParameterProfiles(parameters,objective_function,options)
% [parameters,fh] = getParameterProfiles(...)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least:
%   .number ... number of parameter
%   .guess ... initial guess of parameter
%   .min ... lower bound for parameter values       
%   .max ... upper bound for parameter values       
%   .name = {'name1',...} ... names of the parameters       
%   .MS ... results of global optimization, obtained using for instance 
%       the routine 'getMultiStarts.m'. MS has to contain at least
%       .par ... sorted list n_theta x n_starts of parameter estimates.
%                The first entry is assumed to be the best one.
%       .logPost ... sorted list n_starts x 1 of of log-posterior values
%                corresponding to the parameters listed in .par.
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
%   .parameter_index ... index of the parameters for which the profile
%         is calculated (default = 1:parameters.number).
%   .P.min ... lower bound for profiling parameters, having same
%         dimension as the parameter vector (default = parameters.min).
%   .P.max ... lower bound for profiling parameters, having same
%         dimension as the parameter vector (default = parameters.max).
%   .R_min ... minimal ratio down to which the profile is calculated 
%         (default = 0.03).
%   .dR_max ... maximal relative decrease of ratio allowed
%         for two adjacent points in the profile (default = 0.10) if
%         options.dJ = 0;
%   .dJ ... influnces step size at small likelihood ratio values (default = 0.5).
%   .options_getNextPoint ... options for the generation fo the next profile point
%       .mode ... choice of proposal direction
%           = 'multi-dimensional' (default) ... all parameters are updated.
%               The direct is the same as between the last two points.
%           = 'one-dimensional' ... only parameter for which profile is
%               currently calculated is updated.
%       .guess = 1e-2 ... guess for initial update stepsize
%       .min = 1e-6 ... lower bound for update stepsize
%       .min = 1e2 ... upper bound for update stepsize
%       .update = 1.25 ... incremental change if stepsize is too large or
%           too small.
%   .calc_profiles ... flag for profile calculation
%       = 'true' (default) ... profiles are calculated
%       = 'false' ... profiles are not calculated 
%   .reopt_profil ... flag for profile reoptimization
%       = 'true' ... profiles are reoptimized
%       = 'false' (default) ... profiles are not reoptimized
%   .plot_options ... plot options for plotPropertyProfiles.m.
%   .mode ... output of algorithm
%       = 'visual' (default) ... plots are gnerated which show the progress
%       = 'text' ... optimization results for multi-start is printed on screen
%       = 'silent' ... no output during the multi-start local optimization
%   .fh ... handle of figure in which results are printed. If no
%       handle is provided, a new figure is used.
%   .save ... determine whether results are directly saved
%       = false (default) ... results are not saved
%       = true ... results are stored do an extra folder
%   .foldername ... name of the folder in which results are stored.
%       If no folder is provided, a random foldername is generated.
%   .MAP_index ... index MAP parameter vector starting from which the
%       profile is calculated. This option is helpful if local
%       optima are present.
%
% Outputs:
% ========
% parameters ... updated parameter object containing:
%   .P(i) ... profile for i-th parameter
%       .par ... MAPs along profile
%       .logPost ... maximum log-posterior along profile
%       .R ... ratio
% fh ... figure handle
%
% 2012/05/16 Jan Hasenauer
% 2014/06/12 Jan Hasenauer

% function [parameters,fh] = getParameterProfiles(parameters,objective_function,options)
function [parameters,fh] = getParameterProfiles(varargin)

%% Check and assign inputs
if nargin >= 2
    parameters = varargin{1};
    objective_function = varargin{2};
else
    error('getParameterProfiles requires at least two inputs.')
end

% Check and assign options
options.fmincon = optimset('algorithm','interior-point',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',2000,...
                           'MaxFunEvals',300*parameters.number,...
                           'PrecondBandWidth',inf);
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.mode = 'visual'; % 'text','silent'
options.save = false; % true
options.plot_options.interval = 'dynamic';
options.plot_options.mark_constraint = 'false';
options.plot_options.hold_on = 'false';
options.parameter_index = 1:parameters.number;
options.P.min = parameters.min;
options.P.max = parameters.max;
options.R_min = 0.03;
options.dR_max = 0.1;
options.dJ = 0.5;
options.options_getNextPoint.mode = 'multi-dimensional'; %'one-dimensional';
options.options_getNextPoint.guess = 1e-2;
options.options_getNextPoint.min = 1e-6;
options.options_getNextPoint.max = 1e2;
options.options_getNextPoint.update = 1.25;
options.foldername = strrep(datestr(now,31),' ','__');
options.calc_profiles = 'true'; % 'false'
options.reopt_profiles = 'false'; % 'true'
options.MAP_index = 1;
options.fh = [];
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
        fprintf(' \nProfile likelihood caclulation:\n===============================\n');
    case 'silent' % no output
        % Force fmincon to be silent.
        options.fmincon = optimset(options.fmincon,'display','off');
end

%% Initialization of parameter struct
for i = options.parameter_index
    parameters.P(i).par = parameters.MS.par(:,options.MAP_index);
    parameters.P(i).logPost = parameters.MS.logPost(options.MAP_index);
    parameters.P(i).R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));
end
logPost_max = parameters.MS.logPost(1);

%% Preperation of folder
if options.save
    [~,~,~] = mkdir(options.foldername);
    save([options.foldername '/init'],'parameters');
end

%% Profile calculation -- SEQUENTIAL
if strcmp(options.comp_type,'sequential') && strcmp(options.calc_profiles,'true')

% Profile calculation
for i = options.parameter_index
    % Initialization
    P_par = parameters.MS.par(:,options.MAP_index);
    P_logPost = parameters.MS.logPost(options.MAP_index);
    P_R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));

    % Construction of index set
    I1 = [1:i-1]';
    I2 = [i+1:parameters.number]';
    I  = [I1;I2];

    % Compute profile for in- and decreasing theta_i
    for s = [-1,1]
        % Starting point
        theta  = parameters.MS.par(:,1);
        logPost = parameters.MS.logPost(1);
        
        % Lower and upper bounds for profiles
        theta_min = [parameters.min(I1);options.P.min(i);parameters.min(I2)];
        theta_max = [parameters.max(I1);options.P.max(i);parameters.max(I2)];
        
        % Initialize direction
        dtheta = zeros(parameters.number,1);
        dtheta(i) = s*options.options_getNextPoint.guess;
        
        % Sequential update
        while (options.P.min(i) < theta(i)) && (theta(i) < options.P.max(i)) && ...
              (logPost >= (log(options.R_min) + parameters.MS.logPost(1)))
          
            % Proposal of next profile point
            [theta_next,J_exp] = ...
                getNextPoint(theta,theta_min,theta_max,dtheta/abs(dtheta(i)),...
                             abs(dtheta(i)),options.options_getNextPoint.min,options.options_getNextPoint.max,options.options_getNextPoint.update,...
                             -(log(1-options.dR_max)+options.dJ*(logPost-logPost_max)+logPost),@(theta) obj(theta,objective_function,options.obj_type),...
                             parameters.constraints,options.options_getNextPoint.mode,i);

            % Construction of reduced linear constraints
            [A,b,Aeq,beq] = getConstraints(theta,parameters,I);
            
            % Optimization
            [theta_I_opt,J_opt] = ...
                fmincon(@(theta_I) obj([theta_I(I1);theta_next(i);theta_I(I2-1)],objective_function,options.obj_type,I),... % negative log-posterior function
                                    theta_next(I),...
                                    A  ,b  ,... % linear inequality constraints
                                    Aeq,beq,... % linear equality constraints
                                    parameters.min(I),...   % lower bound
                                    parameters.max(I),...   % upper bound
                                    [],options.fmincon);    % options

            % Restore full vector and determine update direction
            logPost = -J_opt;
            dtheta = [theta_I_opt(I1);theta_next(i);theta_I_opt(I2-1)] - theta;
            theta = theta + dtheta;
            
            % Sorting
            switch s
                case -1
                    P_par = [theta,P_par];
                    P_logPost = [logPost,P_logPost]; 
                    P_R = [exp(logPost - parameters.MS.logPost(1)),P_R];
                case +1
                    P_par = [P_par,theta];
                    P_logPost = [P_logPost,logPost];
                    P_R = [P_R,exp(logPost - parameters.MS.logPost(1))];
            end
            
            % Assignment
            parameters.P(i).par = P_par;
            parameters.P(i).logPost = P_logPost;
            parameters.P(i).R = P_R;

            % Save
            if options.save
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
            end
            
            % Output
            str = [num2str(i,'%d') '-th P: point ' num2str(length(parameters.P(i).R)-1,'%d') ', R = ' ...
                   num2str(exp(- J_opt - parameters.MS.logPost(1)),'%.3e') ' (optimized) / '...
                   num2str(exp(- J_exp - parameters.MS.logPost(1)),'%.3e') ' (predicted)'];
            switch options.mode
                case 'visual', fh = plotParameterProfiles(parameters,'1D',fh,options.parameter_index,options.plot_options); disp(str);
                case 'text', disp(str);
                case 'silent' % no output
            end
        end
    end    
end
end

%% Profile calculation -- PARALLEL
if strcmp(options.comp_type,'parallel') && strcmp(options.calc_profiles,'true')

% Assignement of profile
P = parameters.P;

% Profile calculation
parfor i = options.parameter_index
    % Initialization
    P_par = parameters.MS.par(:,options.MAP_index);
    P_logPost = parameters.MS.logPost(options.MAP_index);
    P_R = exp(parameters.MS.logPost(options.MAP_index)-parameters.MS.logPost(1));

    % Construction of index set
    I1 = [1:i-1]';
    I2 = [i+1:parameters.number]';
    I  = [I1;I2];

    % Compute profile for in- and decreasing theta_i
    for s = [-1,1]
        % Starting point
        theta  = parameters.MS.par(:,options.MAP_index);
        logPost = parameters.MS.logPost(options.MAP_index);
        
        % Lower and upper bounds for profiles
        theta_min = [parameters.min(I1);options.P.min(i);parameters.min(I2)];
        theta_max = [parameters.max(I1);options.P.max(i);parameters.max(I2)];
        
        % Initialize direction
        dtheta = zeros(parameters.number,1);
        dtheta(i) = s*options.options_getNextPoint.guess;
        
        % Sequential update
        while (options.P.min(i) < theta(i)) && (theta(i) < options.P.max(i)) && ...
              (logPost >= (log(options.R_min) + parameters.MS.logPost(1)))
          
            % Proposal of next profile point
            [theta_next,~] = ...
                getNextPoint(theta,theta_min,theta_max,dtheta/abs(dtheta(i)),...
                             abs(dtheta(i)),options.options_getNextPoint.min,options.options_getNextPoint.max,options.options_getNextPoint.update,...
                             -(log(1-options.dR_max)+options.dJ*(logPost-logPost_max)+logPost),@(theta) obj(theta,objective_function,options.obj_type),...
                             parameters.constraints,options.options_getNextPoint.mode,i);

            % Construction of reduced linear constraints
            [A,b,Aeq,beq] = getConstraints(theta,parameters,I);
            
            % Optimization
            [theta_I_opt,J_opt] = ...
                fmincon(@(theta_I) obj([theta_I(I1);theta_next(i);theta_I(I2-1)],objective_function,options.obj_type,I),... % negative log-posterior function
                                    theta_next(I),...
                                    A  ,b  ,... % linear inequality constraints
                                    Aeq,beq,... % linear equality constraints
                                    parameters.min(I),...   % lower bound
                                    parameters.max(I),...   % upper bound
                                    [],options.fmincon);    % options

            % Restore full vector and determine update direction
            logPost = -J_opt;
            dtheta = [theta_I_opt(I1);theta_next(i);theta_I_opt(I2-1)] - theta;
            theta = theta + dtheta;
            
            % Sorting
            switch s
                case -1
                    P_par = [theta,P_par];
                    P_logPost = [logPost,P_logPost]; 
                    P_R = [exp(logPost - parameters.MS.logPost(1)),P_R];
                case +1
                    P_par = [P_par,theta];
                    P_logPost = [P_logPost,logPost];
                    P_R = [P_R,exp(logPost - parameters.MS.logPost(1))];
            end
            
            % Assignment
            P(i).par = P_par;
            P(i).logPost = P_logPost;
            P(i).R = P_R;

            % Save
            if options.save
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__par.csv'],P_par,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__logPost.csv'],P_logPost,'delimiter',',','precision',12);
                dlmwrite([options.foldername '/P' num2str(i,'%d') '__R.csv'],P_R,'delimiter',',','precision',12);
            end
        end
    end    
end

% Assignment
parameters.P = P;

% Output
switch options.mode
    case 'visual', fh = plotParameterProfiles(parameters,'1D',fh,options.parameter_index,options.plot_options);
    case 'text' % no output
    case 'silent' % no output
end

end


% %% REOPTIMIZE PROFILE FROM THE BOARDER
% if strcmp(options.reoptimize,'true')
% disp('');
% disp('Re-optimization of profile');
% % Loop: Parameters
% for i = options.parameter_index
%     disp('');
%     % Index set
%     I1 = [1:i-1]';
%     I2 = [i+1:parameters.number]';
%     I  = [I1;I2];
%     % Likelihood function option
%     options.logPost_options.grad_ind = I(:);
%     options.logPost_options.sign = 'negative';
% 
%     %% COMPUTE OPTIMUM FOR IN-/DECREASING THETA
%     for s = [-1,1]
%         % Find set
%         if s == -1
%             ind = find(parameters.P(i).par(i,:) < parameters.MS.par(i));
%             ind = ind(2:end);
%         else
%             ind = find(parameters.P(i).par(i,:) > parameters.MS.par(i));
%             ind = ind(end-1:-1:1);
%         end
%  
%         % Loop: Index points
%         for k = ind
%             % Starting point
%             theta_next = parameters.P(i).par(:,k+s);
%             theta_next(i) = parameters.P(i).par(i,k);
%             theta_i = theta_next(i);
%             
%             % Modification of options struct
%             options_fmincon = options.fmincon;
%             if isfield(options.fmincon,'TypicalX')
%                 if ~isempty(options.fmincon.TypicalX)
%                     options_fmincon.TypicalX = options_fmincon.TypicalX(I);
%                 end
%             end
%             if isfield(options.fmincon,'FinDiffRelStep')
%                 if ~isempty(options.fmincon.FinDiffRelStep)
%                     options_fmincon.FinDiffRelStep = options_fmincon.FinDiffRelStep(I);
%                 end
%             end
%             
%             % Linear constraints
%             if ~isempty(parameters.constraints.A)
%                 A = parameters.constraints.A(:,I);
%                 b = parameters.constraints.b - 10^-10 - parameters.constraints.A(:,i)*theta_i;
%             else
%                 A = [];
%                 b = [];
%             end
%             if ~isempty(parameters.constraints.Aeq)
%                 Aeq = parameters.constraints.Aeq(:,I) - parameters.constraints.Aeq(:,i)*theta_i;
%                 beq = parameters.constraints.beq;
%             else
%                 Aeq = [];
%                 beq = [];
%             end
%             
%             % Optimize
%             [theta_next,J_opt] = ...
%                 fmincon(@(theta_I) objective_function([theta_I(I1);theta_next(i);theta_I(I2-1)],options.logPost_options),...     % negative log-likelihood function
%                                     theta_next(I),...        % initial parameter
%                                     A  ,b  ,...             % linear inequality constraints
%                                     Aeq,beq,...             % linear equality constraints
%                                     parameters.min(I),...   % lower bound
%                                     parameters.max(I),...   % upper bound
%                                     [],options_fmincon);    % options
%             logPost = -J_opt;
%             
%             % Assignment
%             if -J_opt > parameters.P(i).logPost(k);
%                 parameters.P(i).par(I,k) = theta_next;
%                 parameters.P(i).logPost(k) = logPost;
%                 parameters.P(i).R(k) = exp(logPost - parameters.MS.logPost);
%             end
%             
%             % Update plot
%             if strcmp(options.plot,'true')
%                 fh = plotP(parameters,fh,options.parameter_index,options.plot_options);
%             end
%             
%             % Output command line
%             disp([num2str(i,'%d') '-th P: point ' num2str(k,'%d') ...
%                 ' -> ' num2str(ind(end)-s,'%d')]);
%         end
%     end
%     disp('');
% end
% end
% 
%  

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Profile calculation FINISHED.');
    case 'silent' % no output
end

end


%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
%   I ... index set of optimized parameters
function varargout = obj(theta,fun,type,I)

try
    switch nargout
        case {0,1}
            J = fun(theta);
            if isnan(J)
                error('J is NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J};
                case 'negative log-posterior' , varargout = { J};
            end
        case 2
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J,-G(I)};
                case 'negative log-posterior' , varargout = { J, G(I)};
            end
        case 3
            [J,G,H] = fun(theta);
            if max(isnan([J;G(:);H(:)]))
                error('J, G and/or H contain a NaN.')
            end
            switch type
                case 'log-posterior'          , varargout = {-J,-G(I),-H(I,I)};
                case 'negative log-posterior' , varargout = { J, G(I), H(I,I)};
            end
    end
catch
    switch nargout
        case {0,1}
            varargout = {inf};
        case 2
            varargout = {inf,zeros(length(I),1)};
        case 3
            varargout = {inf,zeros(length(I),1),zeros(length(I))};
    end
end

end

%% Constraint generation
% This function is used to generate the linear constraints for the
% reduced system. 
%   theta ... parameter vector
%   parameter struct ...
%   I ... index set of optimized parameters
function [A,b,Aeq,beq] = getConstraints(theta,parameters,I)

% Index set of parameters which are eliminated
i = setdiff(1:parameters.number,I);

% Reduction of constraints to remaining dimensions
if ~isempty(parameters.constraints.A)
    A = parameters.constraints.A(:,I);
    if ~isempty(i)
        b = parameters.constraints.b - parameters.constraints.A(:,i)*theta(i);
    end
else
    A = [];
    b = [];
end
if ~isempty(parameters.constraints.Aeq)
    Aeq = parameters.constraints.Aeq(:,I);
    if ~isempty(i)
        beq = parameters.constraints.beq - parameters.constraints.Aeq(:,i)*theta(i);
    end
else
    Aeq = [];
    beq = [];
end

end

%% getNextStepProfile is a support function for the profile calculation
%   and is called by computeProfile. It determines the length of the 
%   update step given update direction, parameter constraints,
%   log-posterior and target log posterior.
%
% USAGE:
% ======
% function [theta_next,logPost] = getNextPoint(theta,theta_min,theta_max,dtheta,logPost_target,objective_function)
%
% INPUTS:
% =======
% theta ... starting parameter   
% theta_min ... lower bound for parameters   
% theta_max ... upper bound for parameters   
% dtheta ... upper direction
% logPost_target ... target value for log-posterior
% objective_function ... log-posterior of model as function of the parameters.
%
% Outputs:
% ========
% theta_next ... parameter proposal
% logPost ... log-posterior at proposed parameter
%
% 2012/07/12 Jan Hasenauer

function [theta,J] = getNextPoint(theta,theta_min,theta_max,dtheta,c,c_min,c_max,c_update,J_target,obj,constraints,update_mode,i)

% Initialization
% 1) modification of dtheta
switch update_mode
    case 'multi-dimensional'
        % nothing has to be done
    case 'one-dimensional'
        dtheta([1:i-1,i+1:end]) = 0;
end

% 1) line search
if dtheta(i) > 0 % increasing
    c_bound = (theta_max(i)-theta(i))/dtheta(i);
else
    c_bound = (theta_min(i)-theta(i))/dtheta(i);
end
if c_bound > c_min
    c_max = min(c_max,c_bound);
    c = min(max(c,c_min),c_max);
    search = 1;
else
    c_min = c_bound;
    c_max = c_bound;
    c = c_bound;
    search = 0;
end
% 2) inequality constraints
if ~isempty(constraints.A)
    A = constraints.A;
    b = constraints.b;
else
    A = zeros(1,length(theta));
    b = 1;
end
% 3) parameter projection
theta_fun = @(c) max(min(theta + c*dtheta,theta_max),theta_min);

% Search
theta = theta_fun(c);
if  min(A*theta <= b)
    J = obj(theta);
else
    J = inf;
end

if search == 1
    if J > J_target % => initial c too large
        stop = 0;
        while stop == 0
            c = min(max(c/c_update,c_min),c_max);
            theta = theta_fun(c);
            if c == c_min % lower bound reached
                stop = 1;
                J = obj(theta);
            elseif min(A*theta <= b) % feasible
                J = obj(theta);
                if J <= J_target % objective smaller than target value
                    stop = 1;
                end
            end
        end
    else % => initial c too small
        stop = 0;
        while stop == 0
            cn = min(max(c*c_update,c_min),c_max);
            thetan = theta_fun(cn);
            if min(A*theta <= b) % feasible
                Jn = obj(thetan);
                if Jn <= J_target % objective smaller than target value
                    c = cn;
                    theta = thetan;
                    J = Jn;
                    if cn == c_max % upper bound reached
                        stop = 1;
                    end
                else
                    stop = 1;
                end
            else
                stop = 1;
            end
        end
    end
end
end