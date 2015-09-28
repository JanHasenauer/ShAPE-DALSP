% getParameterSamples.m performs performing a Bayesian uncertainty analysis 
%   using, e.g. adaptive MCMC sampling.
%
% USAGE:
% ======
% [...] = getParameterSamples(parameters,logPosterior)
% [...] = getParameterSamples(parameters,logPosterior,options)
% [parameters] = getParameterSamples(...)
% [parameters,fh_logPost_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis] = getParameterSamples(...)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least:
%   .number ... number of parameter
%   .ml .. maximum likelihood estimate
%   .min ... lower bound for parameter values       
%   .max ... upper bound for parameter values       
% logPosterior ... log-posterior of model as function of the parameters.
% options ... options of algorithm
%   .nsimu_warmup ... length of MCMC warm-up run (default = 1e4).
%   .nsimu_run ... length of MCMC run(default = 5e4).
%   .algorithm ... MCMC sampling scheme (default = 'dram')
%   .qcov ... initial covariance matrix for MCMC sampling
%       (default = 0.001*eye(parameters.number)).
%   .adaptint ... number of function evaluations between adaptation of
%       the MCMC transition kernel (default = 20*parameters.number).
%   .plot ... visualization of the results after the computation (default = 'true').
%   .plot_options ... plot options:
%       .interval ... method uses to determine plot bounds (default = 'dynamic').
%       .hold_on ... conserve of current plot content (default = 'false').
%   .fh_logPost_trace ... figure handle for log-posterior trace plot.
%   .fh_par_trace ... figure handle for parameter trace plots.
%   .fh_par_dis ... figure handle for the parameter distribution plot.
%   .rng ... initialization of random number generator (default = 0).
%       = any ral number r => random generator is initialized with r.
%       = [] ... random number generator is not initialized.
%       (Initializing the random number generator can be helpfult to
%       understand problems.)
%   .plot_options ... plot options for plotPropertyProfiles.m.
%
% Outputs:
% ========
% parameters ... updated parameter object containing:
%   .S ... parameter and posterior sample.
%       .logPost ... log-posterior function along chain
%       .par  ... parameters along chain
% fh_logPost_trace .. figure handle for log-posterior trace
% fh_par_trace .. figure handle for parameter traces
% fh_par_dis .. figure handle for parameter distribution
%
% 2012/07/11 Jan Hasenauer
% 2015/04/29 Jan Hasenauer

% function [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis] = getParameterSamples(parameters,logLikelihood,options)
function [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis_1D,fh_par_dis_2D] = getParameterSamples(varargin)


%% Check and assign inputs
if nargin >= 2
    parameters = varargin{1};
    objective_function = varargin{2};
else
    error('getParameterSamples requires at least three inputs.')
end

% Set options
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.mode = 'visual'; % 'text','silent'
options.save = 'false'; % 'true'
options.rng = 0;
options.thinning = 10;
options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;
options.sampling_scheme = 'DRAM';
options.DRAM.algorithm = 'dram';

options.proposal_scheme = 'MALA'; % 'AM'
options.theta_0 = zeros(parameters.number,1);
options.Sigma_0 = 0.001*eye(parameters.number);

% MALA options
options.MALA.min_regularisation = 1e-6; % minimal regularistion for hessian matrix
options.MALA.w_hist = 0; % interpolation between MALA and AM proposal

% AM options
options.AM.min_regularisation = 1e-6; % minimal regularistion for covariance matrix
options.AM.init_memory_length = 20*parameters.number;
options.AM.adaption_interval = 10*parameters.number;
options.AM.min_acc = 0.15;
options.AM.max_acc = 0.30;
options.AM.adap_Sigma_scale = 0.8;
% Finite memory:
options.AM.adaption_scheme = 'difference';
options.AM.memory_length = 10*parameters.number;

options.SCMC.n_proposals = 10;

% % Classical:
% options.AM.adaption_scheme = 'position';
% options.AM.memory_length = inf;

options.plot_options.interval = 'dynamic';
options.fh.logPost_trace = [];
options.fh.par_trace = [];
options.fh.par_dis_1D = [];
options.fh.par_dis_2D = [];

if isfield(parameters,'MS')
    options.theta_0 = parameters.MS.par(:,1);
    if isfield(parameters.MS,'hessian')
        options.Sigma_0 = inv(parameters.MS.hessian(:,:,1));
    end
else
    if isfield(parameters,'guess')
        options.theta_0 = parameters.guess;
    end
end

if nargin == 3
    options = setdefault(varargin{3},options);
end

%% Initialization of random number generator
if ~isempty(options.rng)
    rng(options.rng);
end

%% Initialization and figure generation
fh_logPost_trace = [];
fh_par_trace = [];
fh_par_dis_1D = [];
switch options.mode
    case 'visual'
        % logL trace
        if isempty(options.fh.logPost_trace)
            fh_logPost_trace = figure;
        else
            fh_logPost_trace = figure(options.fh.logPost_trace);
        end
        % parameter traces
        if isempty(options.fh.par_trace)
            fh_par_trace = figure;
        else
            fh_par_trace = figure(options.fh.par_trace);
        end
        % parameter distribution
        if isempty(options.fh.par_dis_1D)
            fh_par_dis_1D = figure;
        else
            fh_par_dis_1D = figure(options.fh.par_dis_1D);
        end
        if isempty(options.fh.par_dis_2D)
            fh_par_dis_2D = figure;
        else
            fh_par_dis_2D = figure(options.fh.par_dis_2D);
        end
    case 'text'
        fprintf(' \nSampling:\n=========\n');
    case 'silent' % no output
end

%% Selection of sampling proceedure
switch options.sampling_scheme
    
    %% DRAM
    case 'DRAM'
        % This section provides the interface to the MATLAB tooparameters.minox for
        % delayed rejection adaptive metropolis sampling developed by 
        % H. Haario et al. (2006), DRAM: Efficient adaptive MCMC,
        % Stat. Comp., 4(16):339-354.
        
        % Model
        for i = 1:parameters.number
            params{i} = {parameters.name{i},options.theta_0(i),parameters.min(i),parameters.max(i)};
        end

        model.ssfun = @(theta,dummi) 2*logPost(theta,objective_function,options.obj_type,'negative');
        model.sigma2 = 1;
        model.N = 1;

        % Options
        options_dram.method      = options.DRAM.algorithm; % adaptation method (mh,am,dr,dram)
        options_dram.qcov        = options.Sigma_0;      % proposal covariance
        options_dram.adaptint    = options.AM.adaption_interval;  % adaptation interval
        options_dram.printint    = 0;  % how often to show info on acceptance ratios
        options_dram.verbosity   = 0;  % how much to show output in Matlab window
        options_dram.updatesigma = 0;  % update error variance
        options_dram.stats       = 0;  % save extra statistics in results
        if strcmp(options.mode,'visual')
            options_dram.waitbar     = 1;  % show garphical waitbar
        else
            options_dram.waitbar     = 0;
        end

        % Warm-up
        options_dram.nsimu = options.nsimu_warmup; % # simulations
        [results] = mcmcrun(model,[],params,options_dram);

        % Sampling
        options_dram.nsimu = options.nsimu_run; % # simulations
        [~,Theta,~,Obj] = mcmcrun(model,[],params,options_dram,results);

        % Reassignment
        parameters.S.logPost = -0.5*Obj;
        parameters.S.par = Theta';
        
    case 'single-chain'            

        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.nsimu_run),1);
        j = 0;
        acc = 0;
        
        theta = options.theta_0;
        mu_hist = options.theta_0;
        Sigma_hist = options.Sigma_0;
        
        Sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive');
                mu = theta;
                Sigma = options.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive');
                if logP < inf
                [mu,Sigma] = getProposal(theta,G,H,options.MALA.min_regularisation,options.MALA.w_hist,...
                    options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end
        
        % Generate Markov chain
        for i = 1:(options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            if mod(i,100) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                             ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
                    case 'text', disp(str);
                    case 'silent' % no output
                end
            end

            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma)';

            % Evaluate objective function
            if (sum(theta_i < parameters.min) + sum(theta_i > parameters.max) == 0)
                inbounds = 1;
                switch options.proposal_scheme
                    case {'MH','AM'}
                        % Compute log-posterior
                        [logP_i] = logPost(theta_i,objective_function,options.obj_type,'positive');

                        % Update mu and Sigma of proposal
                        mu_i = theta_i;
                        Sigma_i = Sigma;
                    case 'MALA'
                        % Compute log-posterior, gradient and hessian
                        [logP_i,G_i,H_i] = logPost(theta_i,objective_function,options.obj_type,'positive');

                        % Update mu and Sigma of proposal
                        if logP_i < inf
                        [mu_i,Sigma_i] = getProposal(theta_i,G_i,H_i,options.MALA.min_regularisation,options.MALA.w_hist,...
                            mu_hist,Sigma_hist,parameters.min,parameters.max);
                        end
                end
            else
                inbounds = 0;
            end

            % Determine acceptance probability
            if (inbounds == 1) && (logP_i < inf)
                % Transition probabilities
                log_p_forward  = logmvnpdf(theta_i,mu  ,Sigma  );
                log_p_backward = logmvnpdf(theta  ,mu_i,Sigma_i);

                % Acceptance probability
                pacc = exp(logP_i - logP + log_p_backward - log_p_forward);
            else
                pacc = 0;
            end

            % Accept or reject
            r = rand;
            if r <= pacc
                acc    = acc + 1;
                theta  = theta_i;
                dtheta = (theta-mu);
                logP   = logP_i;
                mu     = mu_i;
                Sigma  = Sigma_i; % only for MALA relevant
            else
                dtheta = zeros(size(theta));
            end
            
            % Incremental calculation of mean and covariance
            % (with memory length options.AM.memory_length)
            switch options.AM.adaption_scheme
                case 'position'
                    [mu_hist,Sigma_hist] = updateStatistics(mu_hist,Sigma_hist,theta,max(i,options.AM.init_memory_length),...
                        sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                case 'difference'
                    [Sigma_hist] = updateCovariance(Sigma_hist,dtheta,max(i,options.AM.init_memory_length),...
                        sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    mu_hist = theta;
            end
            
            % Proposal update
            if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
                if acc/i < options.AM.min_acc
                    Sigma_scale = Sigma_scale*options.AM.adap_Sigma_scale;
                elseif acc/i > options.AM.max_acc
                    Sigma_scale = Sigma_scale/options.AM.adap_Sigma_scale;
                end
                Sigma = Sigma_scale*Sigma_hist;
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.AM.min_regularisation*eye(parameters.number);
                end
            end
            
            % Store
            if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                j = j + 1;
                parameters.S.par(:,j) = theta;
                parameters.S.logPost(j) = logP;
            end        
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
    case 'single-chain multi-core'            

        n_proposals = options.SCMC.n_proposals;
        useful_total = 0;
        not_useful_total = 0;
        K_set = [];
       
        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.nsimu_run),1);
        j = 0;
        acc = 0;
        K = 0;

        theta = options.theta_0;
        mu_hist = options.theta_0;
        Sigma_hist = options.Sigma_0;
        
        Sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive');
                mu = theta;
                Sigma = options.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive');
                if logP < inf
                [mu,Sigma] = getProposal(theta,G,H,options.MALA.min_regularisation,options.MALA.w_hist,...
                    options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end
        
        % Generate Markov chain
        i = 1;
        while i <= (options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            if min(mod([i-K:i],100)) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                             ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
                    case 'text', disp(str);
                    case 'silent' % no output
                end
            end

            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma,n_proposals)';

            % Evaluate objective function
            for k = 1:n_proposals
                if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max) == 0)
                    inbounds(k) = 1;
                    switch options.proposal_scheme
                        case {'MH','AM'}
                            % Compute log-posterior
                            [logP_i(k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            mu_i(:,k) = theta_i(:,k);
                            Sigma_i(:,:,k) = Sigma;
                        case 'MALA'
                            % Compute log-posterior, gradient and hessian
                            [logP_i(k),G_i,H_i] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            if logP_i(k) < inf
                                [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),G_i,H_i,options.MALA.min_regularisation,options.MALA.w_hist,mu_hist,Sigma_hist,parameters.min,parameters.max);
                            end
                    end
                else
                    inbounds(k) = 0;
                end
                
                % Determine acceptance probability
                if (inbounds(k) == 1) && (logP_i(k) < inf)
                    % Transition probabilities
                    log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu       ,Sigma         );
                    log_p_backward(k) = logmvnpdf(theta       ,mu_i(:,k),Sigma_i(:,:,k));

                    % Acceptance probability
                    pacc(k) = exp(logP_i(k) - logP + log_p_backward(k) - log_p_forward(k));
                else
                    pacc(k) = 0;
                end

                % Accept or reject
                if rand <= pacc(k)
                    acc_flag(k) = 1;
                else
                    acc_flag(k) = 0;
                end

            end

            % Assignment
            K = min(find(acc_flag));
            if isempty(K)
                K = n_proposals;
            end
            K_set(end+1) = K;
            useful_total = useful_total + K;
            not_useful_total = not_useful_total + n_proposals - K;
            [K,useful_total,not_useful_total]
            for k = 1:K
                % Assignment of new state
                if acc_flag(k) == 1
                    acc    = acc + 1;
                    theta  = theta_i(:,k);
                    dtheta = (theta-mu);
                    logP   = logP_i(k);
                    mu     = mu_i(:,k);
                    Sigma  = Sigma_i(:,:,k); % only for MALA relevant
                else
                    dtheta = zeros(size(theta));
                end
                            
                % Incremental calculation of mean and covariance
                % (with memory length options.AM.memory_length)
                switch options.AM.adaption_scheme
                    case 'position'
                        [mu_hist,Sigma_hist] = updateStatistics(mu_hist,Sigma_hist,theta,max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    case 'difference'
                        [Sigma_hist] = updateCovariance(Sigma_hist,dtheta,max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                        mu_hist = theta;
                end

                % Proposal update
                if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
                    if acc/i < options.AM.min_acc
                        Sigma_scale = Sigma_scale*options.AM.adap_Sigma_scale;
                    elseif acc/i > options.AM.max_acc
                        Sigma_scale = Sigma_scale/options.AM.adap_Sigma_scale;
                    end
                    Sigma = Sigma_scale*Sigma_hist;
                end
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.AM.min_regularisation*eye(parameters.number);
                end
                
                % Store
                if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                    j = j + 1;
                    parameters.S.par(:,j) = theta;
                    parameters.S.logPost(j) = logP;
                end
                
                % Update of counter
                i = i + 1;
            end
            
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
        save('test','K_set');
    case 'multi-chain'

        % Initialization
        temp = linspace(0,1,options.MC.n_temps).^options.MC.exp_temps;

        j = 0;
        parameters.S.PT.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
        parameters.S.PT.logPost = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);

        acc = zeros(options.MC.n_temps,1);
        acc_swap = zeros(options.MC.n_temps-1,1);

        mu = nan(parameters.number,options.MC.n_temps);
        mu_i = nan(parameters.number,options.MC.n_temps);
        mu_hist = repmat(options.theta_0,[1,options.MC.n_temps]);
        Sigma = nan(parameters.number,parameters.number,options.MC.n_temps);
        Sigma_i = nan(parameters.number,parameters.number,options.MC.n_temps);
        Sigma_hist = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);

        Sigma_scale = ones(options.MC.n_temps,1);
        
        logL = nan(options.MC.n_temps,1);
        logL_i = nan(options.MC.n_temps,1);
        logPrior = nan(options.MC.n_temps,1);
        logPrior_i = nan(options.MC.n_temps,1);

        theta = repmat(options.theta_0,[1,options.MC.n_temps]);
        dtheta = zeros(parameters.number,options.MC.n_temps);
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                % Objective function evaluation
                [logL_0,logPrior_0] = ...
                    logPost_pt(options.theta_0,objective_function,options.obj_type,'positive');

                % Assignment
                logL = repmat(logL_0,[options.MC.n_temps,1]);
                logPrior = repmat(logPrior_0,[options.MC.n_temps,1]);
                mu = repmat(options.theta_0,[1,options.MC.n_temps]);
                Sigma = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);
            case 'MALA'
                for k = 1:options.MC.n_temps
                    % Objective function evaluation
                    [logL(k),logPrior(k),dlogL(:,k),dlogPrior(:,k),ddlogL(:,:,k),ddlogPrior(:,:,k)] = ...
                        logPost_pt(options.theta_0,objective_function,options.obj_type,'positive');
                    
                    % Assignment
                    if logL(k) < inf
                    [mu(:,k),Sigma(:,:,k)] = getProposal(options.theta_0,...
                        temp(k)*dlogL(:,k)+dlogPrior(:,k),temp(k)*ddlogL(:,:,k)+ddlogPrior(:,:,k),...
                        options.MALA.min_regularisation,options.MALA.w_hist,...
                        options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                    end
                end
        end
       
        if isnan(logL(end)) || (logL(end) == -inf)
            error(['log-posterior undefined at initial point or error during '...
                   'objective function evaluation in multi-chain sampler. '   ...
                   '(Probably no second input corresponding to temperature.)']);
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end

        % Generate Markov chain
        for i = 1:(options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            if mod(i,100) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                     ' % (acc = ' num2str(100*acc(end)/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
                    case 'text', disp(str);
                    case 'silent' % no output
                end
            end
            
            for k = 1:options.MC.n_temps
                % Propose new parameter vector
                theta_i(:,k) = mvnrnd(mu(:,k),Sigma(:,:,k))';

                % Evaluate objetive function
                if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max)) == 0
                    inbounds = 1;
                    switch options.proposal_scheme
                        case {'MH','AM'}
                            % Objective function evaluation
                            [logL_i(k),logPrior_i(k)] = ...
                                logPost_pt(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            mu_i(:,k) = theta_i(:,k);
                            Sigma_i(:,:,k) = Sigma(:,:,k);
                        case 'MALA'
                            % Objective function evaluation
                            [logL_i(k),logPrior_i(k),dlogL_i(:,k),dlogPrior_i(:,k),ddlogL_i(:,:,k),ddlogPrior_i(:,:,k)] = ...
                                logPost_pt(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Assignment
                            if logL_i(k) < inf
                                logL_i(k)
                                ddlogL_i(:,:,k)
                            [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),...
                                temp(k)*dlogL_i(:,k)+dlogPrior_i(:,k),temp(k)*ddlogL_i(:,:,k)+ddlogPrior_i(:,:,k),...
                                options.MALA.min_regularisation,options.MALA.w_hist,...
                                mu_hist(:,k),Sigma_hist(:,:,k),parameters.min,parameters.max);
                            end
                    end
                else
                    inbounds = 0;
                end
                
                % Determine acceptance probability
                if (inbounds == 1) && (logL_i(k) < inf)
                    % Transition probabilities
                    log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu(:,k)  ,Sigma(:,:,k)  );
                    log_p_backward(k) = logmvnpdf(theta(:,k)  ,mu_i(:,k),Sigma_i(:,:,k));

                    % Acceptance probability
                    pacc(k) = exp(  (temp(k)*logL_i(k)+logPrior_i(k)) ...
                                  - (temp(k)*logL(k)  +logPrior(k)  ) ...
                                  + log_p_backward(k) - log_p_forward(k));
                else
                    pacc(k) = 0;
                end
                
                % Accept or reject
                if rand <= pacc(k)
                    acc(k)       = acc(k) + 1;
                    theta(:,k)   = theta_i(:,k);
                    dtheta(:,k)  = (theta(:,k)-mu(:,k));
                    logL(k)      = logL_i(k);
                    logPrior(k)  = logPrior_i(k);
                    mu(:,k)      = mu_i(:,k);
                    Sigma(:,:,k) = Sigma_i(:,:,k); % only for MALA relevant
                else
                    dtheta(:,k) = 0;
                end
            end
            
            % Loop: Swaps
            for k = 1:options.MC.n_temps-1
%                 pacc_swap(k) = exp(  ...
%                     + (temp(k)*logL(k+1)+logPrior(k+1)) + (temp(k+1)*logL(k  )+logPrior(k  )) ...
%                     - (temp(k)*logL(k  )+logPrior(k  )) - (temp(k+1)*logL(k+1)+logPrior(k+1)));
                pacc_swap(k) = exp((temp(k)-temp(k+1))*(logL(k+1)-logL(k)));

                if rand <= pacc_swap(k)
                    acc_swap(k)        = acc_swap(k) + 1;
                    theta(:,[k,k+1])   = theta(:,[k+1,k]);
                    logL([k,k+1])      = logL([k+1,k]);
                    logPrior([k,k+1])  = logPrior([k+1,k]);
                    mu(:,[k,k+1])      = mu(:,[k+1,k]);
                    switch options.proposal_scheme
                        case 'MALA'
                            Sigma(:,:,[k,k+1]) = Sigma(:,:,[k+1,k]); % only for MALA relevant
                    end
                end
            end
            
            % Loop: Mean, covariance and proposal update
            for k = 1:options.MC.n_temps
                % Incremental calculation of mean and covariance
                % (with memory length options.AM.memory_length)
                switch options.AM.adaption_scheme
                    case 'position'
                        [mu_hist(:,k),Sigma_hist(:,:,k)] = updateStatistics(mu_hist(:,k),Sigma_hist(:,:,k),theta(:,k),max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    case 'difference'
                        [Sigma_hist(:,:,k)] = updateCovariance(Sigma_hist(:,:,k),dtheta(:,k),max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                        mu_hist(:,k) = theta(:,k);
                end
            
                % Proposal update
                if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
%                     Sigma_scale(k) = Sigma_scale(k) * (acc(k)/i)/(options.AM.min_acc + options.AM.max_acc);
                    if acc(k)/i < options.AM.min_acc
                        Sigma_scale(k) = Sigma_scale(k)*options.AM.adap_Sigma_scale;
                    elseif acc(k)/i > options.AM.max_acc
                        Sigma_scale(k) = Sigma_scale(k)/options.AM.adap_Sigma_scale;
                    end
                    Sigma(:,:,k) = Sigma_scale(k)*Sigma_hist(:,:,k);

                    % Regularisation
                    [~,p] = cholcov(Sigma(:,:,k),0);
                    if p ~= 0
                        Sigma(:,:,k) = Sigma(:,:,k) + options.AM.min_regularisation*eye(parameters.number);
                    end
                end
            end
            
            % Store
            if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                j = j + 1;
                parameters.S.PT.par(:,j,:) = theta;
                parameters.S.PT.logPost(j,:) = logL + logPrior;
                
                str = num2str(100*acc(1)/i,'%.2f');
                for k = 2:options.MC.n_temps
                    str = [str ', ' num2str(100*acc(k)/i,'%.2f')];
                end
                disp(['acc = ' str ' %']);

                str = num2str(100*acc_swap(1)/i,'%.2f');
                for k = 2:options.MC.n_temps-1
                    str = [str ', ' num2str(100*acc_swap(k)/i,'%.2f')];
                end
                disp(['acc_swap = ' str ' %']);
            end
        end
        
        % Reduction
        parameters.S.PT.par = parameters.S.PT.par(:,1:j,:);
        parameters.S.PT.logPost = parameters.S.PT.logPost(1:j,:);
        parameters.S.par = parameters.S.PT.par(:,1:j,end);
        parameters.S.logPost = parameters.S.PT.logPost(1:j,end);

        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end

%         % Check
%         try
%             objective_function(options.theta_0);
%         catch
%             error(['Error during objective function evaluation in multi-chain sampler. '  ...
%                    '(Probably no second input corresponding to temperature.)']);
%         end
%             
%         % Initialization
%         temp = linspace(0,1,options.MC.n_temps).^options.MC.exp_temps;
% 
%         parameters.S.PT.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
%         parameters.S.PT.logPost = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
%         j = 0;
%         acc = zeros(options.MC.n_temps,1);
%         acc_trans = zeros(options.MC.n_temps-1,1);
% 
%         theta = repmat(options.theta_0,[1,options.MC.n_temps]);
%         mu_hist = repmat(options.theta_0,[1,options.MC.n_temps]);
%         Sigma_hist = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);
%         
%         Sigma_scale = ones(options.MC.n_temps,1);
%         
%         logP = nan(options.MC.n_temps,1);
%         logP_i = nan(options.MC.n_temps,1);
%         dtheta = nan(parameters.number,options.MC.n_temps);
%         
%         % Initialization and testing of starting point
%         for k = 1:options.MC.n_temps
%             switch options.proposal_scheme
%                 case {'MH','AM'}
%                     [logP(k)] = logPost(theta(:,k),objective_function,options.obj_type,'positive',temp(k));
%                     mu(:,k) = theta(:,k);
%                     Sigma(:,:,k) = options.Sigma_0;
%                 case 'MALA'
%                     [logP(k),G(:,k),H(:,:,k)] = logPost(theta(:,k),objective_function,options.obj_type,'positive',temp(k));
%                     [mu(:,k),Sigma(:,:,k)] = getProposal(theta(:,k),G(:,k),H(:,:,k),options.MALA.min_regularisation,options.MALA.w_hist,...
%                         options.theta_0,options.Sigma_0,parameters.min,parameters.max);
%             end
%         end
%         if isnan(logP(end)) || (logP(end) == -inf)
%             error('log-posterior undefined at initial point.');
%         end
% 
%         % Initialization of waitbar
%         h = waitbar(0,['MCMC sampling completed to 0 % (acc = 0 %)']);
% 
%         % Generate Markov chain
%         for i = 1:(options.nsimu_run+options.nsimu_warmup)
%             % Report of progress
%             if mod(i,100) == 0
%                 waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,...
%                     ['MCMC sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
%                      ' % (acc = ' num2str(100*acc(end)/i,'%.2f') ' % )']);
%             end
% 
%             for k = 1:options.MC.n_temps
%                 % Propose new parameter vector
%                 theta_i(:,k) = mvnrnd(mu(:,k),Sigma(:,:,k))';
% 
%                 % Check bounds
%                 if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max) == 0)
%                     switch options.proposal_scheme
%                         case {'MH','AM'}
%                             % Compute log-posterior
%                             [logP_i(k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',temp(k));
% 
%                             % Update mu and Sigma of proposal
%                             mu_i(:,k) = theta_i(:,k);
%                             Sigma_i(:,:,k) = Sigma(:,:,k);
%                         case 'MALA'
%                             % Compute log-posterior, gradient and hessian
%                             [logP_i(k),G_i(:,k),H_i(:,:,k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',temp(k));
% 
%                             % Update mu and Sigma of proposal
%                             [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),G_i(:,k),H_i(:,:,k),options.MALA.min_regularisation,options.MALA.w_hist,...
%                                 mu_hist(:,k),Sigma_hist(:,:,k),parameters.min,parameters.max);
%                     end
% 
%                     % Transition probabilities
%                     log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu(:,k)  ,Sigma(:,:,k)  );
%                     log_p_backward(k) = logmvnpdf(theta(:,k)  ,mu_i(:,k),Sigma_i(:,:,k));
% 
%                     % Acceptance probability
%                     pacc(k) = exp(logP_i(k) - logP(k) + log_p_backward(k) - log_p_forward(k));
%                 else
%                     pacc(k) = 0;
%                 end
% 
%                 % Accept or reject
%                 r(k) = rand;
%                 if r(k) <= pacc(k)
%                     acc(k)       = acc(k) + 1;
%                     theta(:,k)   = theta_i(:,k);
%                     dtheta(:,k)  = (theta(:,k)-mu(:,k));
%                     logP(k)      = logP_i(k);
%                     mu(:,k)      = mu_i(:,k);
%                     Sigma(:,:,k) = Sigma_i(:,:,k); % only for MALA relevant
%                 else
%                     dtheta(:,k) = 0;
%                 end
%                 
%                 % Incremental calculation of mean and covariance
%                 % (with memory length options.AM.memory_length)
%                 switch options.AM.adaption_scheme
%                     case 'position'
%                         [mu_hist(:,k),Sigma_hist(:,:,k)] = updateStatistics(mu_hist(:,k),Sigma_hist(:,:,k),theta(:,k),max(i,options.AM.init_memory_length),...
%                             sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
%                     case 'difference'
%                         [Sigma_hist(:,:,k)] = updateCovariance(Sigma_hist(:,:,k),dtheta(:,k),max(i,options.AM.init_memory_length),...
%                             sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
%                         mu_hist(:,k) = theta(:,k);
%                 end
%             
%                 % Proposal update
%                 if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
%                     if acc/i < options.AM.min_acc
%                         Sigma_scale(k) = Sigma_scale(k)*options.AM.adap_Sigma_scale;
%                     elseif acc/i > options.AM.max_acc
%                         Sigma_scale(k) = Sigma_scale(k)/options.AM.adap_Sigma_scale;
%                     end
%                     Sigma(:,:,k) = Sigma_scale(k)*Sigma_hist(:,:,k);
% 
%                     % Regularisation
%                     [~,p] = cholcov(Sigma(:,:,k),0);
%                     if p ~= 0
%                         Sigma(:,:,k) = Sigma(:,:,k) + options.AM.min_regularisation*eye(parameters.number);
%                     end
%                 end
%             end
%             
%             % Store
%             if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
%                 j = j + 1;
%                 parameters.S.PT.par(:,j,:) = theta;
%                 parameters.S.PT.logPost(j,:) = logP;
%                 
%                 str = num2str(100*acc(1)/i,'%.2f');
%                 for k = 2:options.MC.n_temps
%                     str = [str ', ' num2str(100*acc(k)/i,'%.2f')];
%                 end
%                 disp(['acc = ' str ' %']);
% 
%             end        
%         end
%         % Reduction
%         parameters.S.PT.par = parameters.S.PT.par(:,1:j,:);
%         parameters.S.PT.logPost = parameters.S.PT.logPost(1:j,:);
%         parameters.S.par = parameters.S.PT.par(:,1:j,end);
%         parameters.S.logPost = parameters.S.PT.logPost(1:j,end);
end

%% Visualization of results
if strcmp(options.mode,'visual')
    try
        % Diagnosis plots
        plotMCMCdiagnosis(parameters,'log-posterior',fh_logPost_trace);
        plotMCMCdiagnosis(parameters,'parameters',fh_par_trace);

        % Parameter distribution
        plotParameterSamples(parameters,'1D',fh_par_dis_1D,[],options.plot_options);
        %plotParameterSamples(parameters,'2D',fh_par_dis_2D,[],options.plot_options);

        % Chain statistics
        chainstats(parameters.S.par');
    end
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Sampling FINISHED.');
    case 'silent' % no output
end

end



%% Proposal calculating function
% This function determines the mean and covariance of the MALA proposal and
% regularised / adaptive variants of it.
%   theta ... current parameter vector
%   grad ... gradient of log-posterior
%   H ... hessian or hessian approximation of log-posterior
%   beta ... regularisation parameter
%   w_hist ... weighting of history
%       = 0 ... only MALA
%       = 1 ... only adaptive Metropolis
%   mu_hist ... temporal average of theta
%   Sigma_hist ... temporal covariance of theta
%   lb ... lower bound for theta
%   ub ... upper bound for theta
function [mu,Sigma] = getProposal(theta,grad,H,beta,w_hist,mu_hist,Sigma_hist,lb,ub)

% Planning of MALA step
if w_hist ~= 1
    % Regularisation
    [~,p] = cholcov(-H,0);
    if p ~= 0
        k = 0;
        while p ~= 0
            H_k = H - 10^k*beta*eye(length(theta));
            [~,p] = cholcov(-H_k,0);
            k = k+1;
        end
        H = H_k;
    end
    
    % Newton step
    Sigma_MALA = -inv(H);
    Sigma_MALA = 0.5*(Sigma_MALA+Sigma_MALA');
    mu_MALA = theta - H\grad;
end

% Interpolation between
% a)  MALA (w_hist = 0) and
% b)  adaptive Metropolis with stabilized mean (w_hist = 1) 
if w_hist == 0
    % a) MALA
    Sigma = Sigma_MALA;
    mu = mu_MALA;
elseif w_hist == 1
    % b) AM
    Sigma = Sigma_hist;
    mu = mu_hist;
else
    % c) Hybrid of MALA and AM
    Sigma = inv((1-w_hist)*inv(Sigma_MALA) + w_hist*inv(Sigma_hist));
    Sigma = 0.5*(Sigma+Sigma');
    mu = Sigma*((1-w_hist)*inv(Sigma_MALA)*mu_MALA + w_hist*inv(Sigma_hist)*mu_hist);
end

end


%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = logPost(theta,fun,type,sign)

switch sign
    case 'negative'
        s = -1;
    case 'positive'
        s = +1;
end

try
    switch nargout
        case 1
            J = fun(theta);
            if isnan(J)
                error('J is NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J};
                case 'negative log-posterior' , varargout = {s*-J};
            end
        case 2
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:)};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:)};
            end
        case 3
            [J,G,H] = fun(theta);
            if max(isnan([J;G(:);H(:)]))
                error('J, G and/or H contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:),s* H};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:),s*-H};
            end
    end
catch error_msg
    disp(['Objective function evaluation failed because: ' error_msg.message]);
    switch nargout
        case 1
            varargout = {-s*inf};
        case 2
            varargout = {-s*inf,zeros(length(theta),1)};
        case 3
            varargout = {-s*inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Objetive function interface for parallel tempering
% This function is used as interface to the user-provided objective
% function for the parallel tempering sampling algorithm. It adapts the 
% sign and supplies the correct number of outputs. Furthermore, it catches
% errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = logPost_pt(theta,fun,type,sign)

switch sign
    case 'negative'
        s = -1;
    case 'positive'
        s = +1;
end

try
    switch nargout
        case 2
            [logL,logPrior] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {s* logL,s* logPrior};
                case 'negative log-posterior' , varargout = {s*-logL,s*-logPrior};
            end
        case 6
            [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {s* logL,s* logPrior,s* dlogL,s* dlogPrior,s* ddlogL,s* ddlogPrior};
                case 'negative log-posterior' , varargout = {s*-logL,s*-logPrior,s*-dlogL,s*-dlogPrior,s*-ddlogL,s*-ddlogPrior};
            end
    end
catch error_msg
    disp(['Objective function evaluation failed because: ' error_msg.message]);
    switch nargout
        case 2
            varargout = {-s*inf,-s*inf};
        case 6
            varargout = {-s*inf,-s*inf,zeros(length(theta),1),zeros(length(theta),1),zeros(length(theta)),zeros(length(theta))};
    end
end

end
