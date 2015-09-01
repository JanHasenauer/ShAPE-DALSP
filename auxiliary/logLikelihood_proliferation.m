% [log,grad] = logLikelihood_proliferation(theta,M,D,options)
function varargout = logLikelihood_proliferation(varargin)

%% CHECK/ASSIGN INPUTS:
if nargin >= 3
	theta = varargin{1};
	M = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

% Assign defaults for options
options.costfuntype = 'least-squares';% 'approx logL';
options.simulation.noise.flag = 'no';
options.simulation.plot = 'off';
options.fh = [];
options.sign = 'positive';
options.grad_ind = [1:length(theta)]';
options.simulation.dependency = M.type; 
if nargin == 4
	options = setdefault(varargin{4},options);
end

% Get global timer for visualization
if ~isempty(options.fh)
    global tp 
end

%% SIMULATION OPTIONS
% Options
if nargout == 2
    options.simulation.gradient = 'yes';
else
    options.simulation.gradient = 'no';
end

%% SIMULATION OPTIONS
x = [D(1).bins(1,1); D(1).bins(:,2)]';
[Sim] = CPsimulateDALSP(M,theta,D.t,options.simulation.a_sim,x,options.simulation,'estimation');
            
%% COMPUTATION OF BIN-PROBABILITY
H = bsxfun(@times,Sim.n_y(:,1:end-1) + Sim.n_y(:,2:end), 0.5*diff(x(:)'));
if nargout == 2
    dHdtheta = squeeze(bsxfun(@times,Sim.dn_ydtheta(:,1:end-1,:) + Sim.dn_ydtheta(:,2:end,:), 0.5*diff(x(:)')));
end

%% EVALUATION OF LIKELIHOOD
switch options.costfuntype
    % Exact likelihood using multinomial distribution
    case 'exact logL'
        error('This option is not available so far.');
        
    % least-squares
    case 'least-squares'
        % Initialization
        logL = 0;
        if nargout == 2
            grad = zeros(1,length(options.grad_ind));
        end
        % Loop: points in time
        for k = 1:length(D.t)
            logL  =  logL - sum((D.cellcount(k,:)'-H(k,:)').^2,1);
            if nargout == 2
                grad = grad + 2*sum(bsxfun(@times,D.cellcount(k,:)'-H(k,:)',squeeze(dHdtheta(k,:,options.grad_ind))),1);
            end
        end

    % Approximate likelihood using multi-variate Gaussian distribution
    case 'approx logL'
        % Initialization
        logL = 0;
        % Loop: points in time
        for k = 1:length(D.t)
            % Idea of the approximation:
            % Instead of using a binomial model a 
            
            % Sums for normalization
            sumd = sum(D.cellcount(k,:)');
            sump = sum(H(k,:)');
            % Probability distribution
            pn = H(k,:)'/sump;
            % Addition of background
            background = 0.01;
            pn = (pn + background/length(pn))/(1+background);
            % Evaluation of approximated likelihood function:
            %   Sample distribution
            logL = logL - 0.5*sum(log(2*pi*pn*sumd)) ...
                  - 0.5*norm((D.cellcount(k,:)' - sumd*pn)./sqrt(pn*sumd),2)^2;
            %   Size of population
            logL = logL - 0.5*log(2*pi*(options.noise.sigma_mean*sump)^2) ...
                  - 0.5*((sumd-sump)/(options.noise.sigma_mean*sump))^2;
        end
        
    case 'logL'
        % Initialization
        sigma = M.noise_N.sigma_fun(theta);
        dsigmadtheta = M.noise_N.dsigmadtheta_fun(theta);
        dsigmadtheta = dsigmadtheta(:)';
        logL = 0;
        if nargout == 2
            grad = zeros(1,length(options.grad_ind));
        end
        
        % Loop: points in time
        for k = 1:length(D.t)
            % Sums for normalization
            d = D.cellcount(k,:);
            Sd = sum(d);
            Sp = sum(H(k,:));
            if nargout == 2
                dSpdtheta = sum(squeeze(dHdtheta(k,:,options.grad_ind)),1);
            end
            % Probability distribution
            background = 0.01;
            p = (H(k,:)'/Sp + background/length(H(k,:)'))/(1+background);
            if nargout == 2
                dpdtheta = 1/(1+background)*(...
                      squeeze(dHdtheta(k,:,options.grad_ind))/Sp ...
                    - H(k,:)'/(Sp^2)*dSpdtheta ...
                        );
            end
            % Evaluation of likelihood function:
            %   Sample distribution
            I = log(1:Sd);
            logL = logL + sum(I);
            for i = 1:length(d)
                logL = logL - sum(I(1:d(i)));
            end
            logL = logL + d*log(p);
            if nargout == 2
                grad = grad + (D.cellcount(k,:)./p')*dpdtheta;
            end
            %   Size of population
            logL = logL - 0.5*log(2*pi) - log(sigma) - log(Sd) ...
                        - 0.5*((log(Sd)-log(Sp))^2)/(sigma^2);
            if nargout == 2
                grad = grad - dsigmadtheta(options.grad_ind)/sigma ...
                            + (log(Sd)-log(Sp))/(sigma^2)*(dSpdtheta/Sp) ...
                            + ((log(Sd)-log(Sp))^2)/(sigma^3)*dsigmadtheta(options.grad_ind);
            end
        end
       
        
    otherwise
        error('This option is not available.');
end

%% PLOT CURRENT FIT
if ~isempty(options.fh)
    if etime(clock,tp) > 10
        % Reset tp
        tp = clock;
        % Plot
        plotProliferationAssay(Sim,D,options.fh);
    end
end

%% ASSIGN OUTPUT
switch nargout
    % One output
    case {0,1}
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
            case 'negative'
                varargout{1} = -logL;
        end
    % Two outputs
    case 2
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
                varargout{2} =  grad';
            case 'negative'
                varargout{1} = -logL;
                varargout{2} = -grad';
        end
end

% END
end
