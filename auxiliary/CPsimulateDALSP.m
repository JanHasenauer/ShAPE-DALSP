% CPsimulateDALSP.m simulates a division-, age- and label-structured population model 
% using the coupled system of partial differential equations:
%
% \frac{\partial N_0(t,a,x)}{\partial t} + \frac{\partial N_0(t,a,x)}{\partial a} + 
%               + \frac{\partial (\nu(x) N_0(t,a,x))}{\partial x}
%               = - \left( \alpha_0(t,a) + \beta_0(t,a) \right) N_0(t,a,x)
%
%                                       IC: N_0(0,a,x) \equiv N_{0,init}(a,x),
%                                       BC: N_0(t,0,x) \equiv 0
% \frac{\partial N_i(t,a,x)}{\partial t} + \frac{\partial (N_i(t,a,x))}{\partial a}
%               + \frac{\partial (\nu(x) N_i(t,a,x))}{\partial x}
%               = - \left( \alpha_i(t,a) + \beta_i(t,a) \right) N_i(t,a,x) 
%
%                                       IC: N_i(0,a,x) \equiv N_{i,init}(a,x)
%                                                      \equiv 0,
%                                       BC: N_i(t,0,x) \equiv 2 \gamma \cdot
%                                               \int_{\Rm_+} \alpha{i-1}(t,a) \cdot
%                                               N{i-1}(t,a, \gamma x) da ;
%                                        i = 1, \ldots, S
%
% in which N_i(t,a,x) is the number density of label 'x' and cell age 'a' within 
% generation i (this means that the cells divided i times) at time 't'.
%
% USAGE:
% =======
% [ ] = CPsimulateDALSP(M,theta,t,a,x,options)
%
% INPUT:
% =====
% M ... model specification:
%   .alpha ... vector of rates of cell division. The ith entry of vector  
%              alpha provides the rate at which generation i-1 divides. 
%              (The original population is generation 0 and is assumed to 
%              have divided 0 times.) If alpha is scalar, the division rate 
%              is assumed to be constant over the generations.
%   .beta ...  vector of rates of cell death. The ith entry of vector beta 
%			   provides the rate at which generation i-1 undergoes cell 
%			   death. If beta is scalar, the rates of cell death are 
%			   assumed to be constant over the generations.
%   .gamma ... label dilution factor upon cell division.
%   .k ... scalar providing the rate of label degradation.
%   .c ... scalar providing the rate of decay of the label degradation.
%
% t ... vector of points in time at which the distribution within the 
%   population is evaluated.
% a ... cell-ages at which the distribution of the subpopulations is evaluated.
% x ... label concentrations at which the distribution of the subpopulations  
%   is evaluated.
% options ... options of the algorithm:
%	.S ... number of generations simulated (default = 10).
%   .noise ... properties of measurement noise
%      .flag = {'yes','no' (default)} ... flag indicating whether
%               noise is added.
%      .mu ... mu of log-normally distributed measurement noise.
%      .sigma ... sigma of log-normally distributed measurement noise.
%      .shift ... vector multiplicative shift of whole histrogram.
%               (This can occur due to the addition of different amounts of
%               CFSE in the beginning of the experiment. shift = [1,...,1]
%               means that the added amount of label was always identical,
%               while for shift = [1,1.2,1,...,1], the amount of label added
%               to generate the second measurement was 20 % higher.)
%
% OUTPUT:
% =======
% DALSP ... simulation results for division- and label-structured population
%		.alpha ... rates of cell division
%		.beta ... rates of cell death
%		.gamma ... label dilution factor
%		.k ... rate of label dilution
%		.c ... rate of reduction of label dilution
%		.S ... number of subpopulations
%		.t ... t
%		.N ... matrix containing the unnormalized label distributions 
%           within the subpopulations. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different cell ages a
%               3rd dimension: different label concentrations x  
%               4th dimension: different subpopulations i
%		.Nx ... matrix containing the age distributions
%           within the subpopulations. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different cell ages a
%               3rd dimension: different subpopulations i 
%		.Na ... matrix containing the label distributions
%           within the subpopulations. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x 
%               3rd dimension: different subpopulations i 
%		.Nay ... matrix containing the label distributions
%           within the subpopulations with noise added. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x 
%               3rd dimension: different subpopulations i 
%		.Nxa ... matrix containing the number of cells
%           within the subpopulations. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different subpopulations i 
%		.na ... matrix containing the normalized label distributions 
%           within the subpopulations. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x 
%               3rd dimension: different subpopulations i
%		.nay ... matrix containing the normalized label distributions 
%           within the subpopulations with noise added. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x + noise
%               3rd dimension: different subpopulations i
%		.M ... matrix containing the unnormalized age-label distributions 
%           of the overall cell population. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different cell ages a
%               3rd dimension: different label concentration x
%		.Mx ... matrix containing the unnormalized age distributions 
%           of the overall cell population. 
%           Structure:
%               1st dimension: different time points t
%               2nd dimension: different cell ages a
%		.Ma ... matrix containing the unnormalized label distributions 
%           of the overall cell population. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x
%       .May ... matrix containing the unnormalized label distributions 
%           of the overall cell population with noise added. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x
%		.Mxa ... matrix containing the number of cells 
%           of the overall cell population. 
%           Structure:
%               1st dimension: different time points t 
%		.ma ... matrix containing the normalized label distributions 
%           of the overall cell population. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x
%		.may ... matrix containing the normalized label distributions 
%           of the overall cell population with noise added. 
%           Structure:
%               1st dimension: different time points t 
%               2nd dimension: different label concentration x
%       .avediv ... average division number within the cell population.
%           Structure:
%               1st dimension: different time points t
%       .x ... fluorescence value at which the concentration is evaluated.
%
% function [DALSP] = 
%   CPsimulateDALSP(M,theta,t,a,x,options,type)
function [Sim,output,doutput] = CPsimulateDALSP(varargin)

%% CHECK/ASSIGN INPUTS:
if nargin >= 5
    M = varargin{1};
    theta = varargin{2};
	t = varargin{3};  t = t(:);
    a = varargin{4};  a = a(:);
	x = varargin{5};  x = x(:);
else
	error('Not enough inputs!');
end

% Type of use ('simulation' or 'estimation')
if nargin == 7
    type = varargin{7};
else
    type = 'simulation';
end

% Assign options if necessary
switch type
    case 'simulation'
        % Default options
        options.gradient = 'no';
        options.noise.flag = 'no';
        options.plot = 'off';
        if isfield(M,'type')
            options.dependency = M.type; 
        else
            options.dependency = 'time- and age-dependent'; 
        end
        % Check which options are set
        if nargin >= 6
            options = setdefault(varargin{6},options);
        end
        % Check for required variables
        if ~isfield(options,'t_sim') || ...
           ~isfield(options,'a_sim')
            error('Simulation requires at least options: t_sim and a_sim.');
        end
        % Check that time and age increment is identical and constant
        D = [diff(options.t_sim(:));diff(options.a_sim(:))];
        if max(D) > 1.0001*min(D) % Complex check to ensure numerical robustness
            error(['Time and age increments in simulation vectors have to ' ...
                   'be identical and constant.']);
        end
    case 'estimation'
        % It is assumed that all necessary information are provided!
        % Concistency of the inputs in not check to reduce computational
        % complextity.
        options = varargin{6};
    otherwise
        error('This tye of use is unknown.');
end

% Assign variables
S = M.S;
t_sim = options.t_sim(:);
a_sim = options.a_sim(:);
dt = t_sim(2)-t_sim(1);

% Assignment of variable dimensions
dim_x = length(x);
dim_t = length(t);
dim_a = length(a);
dim_theta = length(theta);
dim_t_sim = length(t_sim);
dim_a_sim = length(a_sim);

%% SIMULATE CELL POPULATION
% Due to the properties of the system of PDES described above, the solution of
% the individual subpopulations can be computed extremely efficiently.
% The computation approach we developed is based on the decomposition
% of the solution of the original system of PDEs into a model for the age 
% distribution of cells 'N^x_i(t,a)', and the normalized distribution of label
% 'n^a_i(t,x)' within the generations.
% The complete solution is then derived by
%	'N_i(t,a,x) = N^x_i(t,a) n_i(t,x), \quad i = 1,\ldots S.
% The solution 'N^x_i(t,a)' follows a system of coupled partial differential
% equations, whereas 'n^a_i(t,x)' solves a single partial differential equation
% (PDE). Both solutions can always be computed analytically using the method of
% characteristics.
% (Remarks: N^x_i(t,a) is the integral of N_i(t,a,x) over x.
%           The PDE for N^x_i(t,a) equals the 'von Foerster equation'.)

%% INITIALIZATION
if ~strcmp(type,'estimation')
    Sim.t = t;
    Sim.a = a;
    Sim.x = x;
    Sim.n_ai = zeros(dim_t,dim_a,S);
    Sim.p_xi = zeros(dim_t,dim_x,S);
    if strcmp(options.gradient,'yes')
        Sim.dn_aidtheta = zeros(dim_t,dim_a,S,dim_theta);
        Sim.dp_xidtheta = zeros(dim_t,dim_x,S,dim_theta);
    end
    Sim.b_i  = zeros(dim_t,S);
end

Sim.N_i  = zeros(dim_t,S);
Sim.p_yi = zeros(dim_t,dim_x,S);
if strcmp(options.gradient,'yes')
    Sim.dN_idtheta  = zeros(dim_t,S,dim_theta);
    Sim.dp_yidtheta = zeros(dim_t,dim_x,S,dim_theta);
    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
end

%% MATRICES FOR INTERPOLATION AND INTEGRATION
% Interpolation matrix for time vector
Tr_t = sparse(zeros(dim_t,dim_t_sim));
Tr_ia0 = sparse(zeros(dim_t,dim_a_sim));
Tr_ia  = sparse(zeros(dim_t,dim_a_sim));
for k = 1:dim_t
    % Determine index and weight for interpolation
    ind(1) = find(t_sim <= t(k),1,'last');
    if t_sim(ind(1)) == t(k)
        ind(2) = ind(1);
        Tr_t(k,ind(1)) = 1;
    else
        ind(2) = ind(1)+1;
        Tr_t(k,ind(1)) = (t_sim(ind(2)) - t(k))/diff(t_sim(ind));
        Tr_t(k,ind(2)) = (t(k) - t_sim(ind(1)))/diff(t_sim(ind));
    end
    % Determine weight for integration
    if dim_a_sim-ind(2) >= 2
        Tr_ia0(k,ind(2):dim_a_sim) = dt*[0.5,ones(1,dim_a_sim-ind(2)-1),0.5];
    end    
    if ind(2) >= 2
        Tr_ia(k,1:ind(2)) = dt*[0.5,ones(1,ind(2)-2),0.5];
    end    
end
% Integration matrix for simulation
Tr_ia_sim0 = sparse(triu(dt*ones(dim_t_sim,dim_a_sim),1));
Tr_ia_sim0(1:dim_t_sim+1:min(dim_t_sim^2,dim_t_sim*dim_a_sim)) = 0.5*dt;
Tr_ia_sim0(:,dim_a_sim) = [0.5*dt*ones(min(dim_t_sim,dim_a_sim-1),1);...
                           zeros(dim_t_sim-min(dim_t_sim,dim_a_sim-1),1)];
Tr_ia_sim = sparse(tril(dt*ones(dim_t_sim,dim_a_sim),-1));
Tr_ia_sim(1:dim_t_sim+1:min(dim_t_sim^2,dim_t_sim*dim_a_sim)) = 0.5*dt;
Tr_ia_sim(:,1) = [0;0.5*dt*ones(dim_t_sim-1,1)];
Tr_ia_sim(:,dim_a_sim) = [zeros(min(dim_t_sim,dim_a_sim-1),1);...
                          0.5*dt*ones(dim_t_sim-min(dim_t_sim,dim_a_sim-1),1)];

% Interpolation matrix for age vector
if ~strcmp(type,'estimation')
    Tr_a = sparse(zeros(dim_a_sim,dim_a));
    for k = 1:dim_a
        % Determine index and weight for interpolation
        ind(1) = find(a_sim <= a(k),1,'last');
        ind(2) = min(ind(1)+1,dim_a_sim);
        % ind = [Ind(k),min(Ind(k)+1,dim_a_sim)];
        if a_sim(ind(1)) == a(k)
            Tr_a(ind(1),k) = 1;
        else
            ind(2) = ind(1)+1;
            Tr_a(ind(1),k) = (a_sim(ind(2)) - a(k))/diff(a_sim(ind));
            Tr_a(ind(2),k) = (a(k) - a_sim(ind(1)))/diff(a_sim(ind));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLUTION OF DASP MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISTINCTION BETWEEN CASE WITH 
% a) initial age distribution is available
% b) initial age distribution is not available
%       => all cells start with age zero

if strcmp(M.IC.na0.type,'logn') % a) Solution of DASP model for i = 0
    %% Computation of exponential
    switch options.dependency
        % Time- and age-dependent division and/or death rates
        case 'time- and age-dependent'
            % Definition of grid
            T = repmat(t_sim ,1,dim_a_sim);
            A = repmat(a_sim',dim_t_sim,1);
            % Evaluation of alpha and beta
            alpha_i = reshape(M.alpha_fun{1}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
            beta_i  = reshape( M.beta_fun{1}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
            if strcmp(options.gradient,'yes')
                dalpha_idtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                dbeta_idtheta  = zeros(dim_t_sim,dim_a_sim,dim_theta);
                for l = 1:dim_theta
                    dalpha_idtheta(:,:,l) = reshape(M.dalphadtheta_fun{1,l}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
                    dbeta_idtheta(:,:,l)  = reshape(M.dbetadtheta_fun{1,l}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
                end
            end
            % Computation of exponential
            expI = zeros(dim_t_sim,dim_a_sim);
            if strcmp(options.gradient,'yes')
                dexpIdtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
            end
            for j = 1:dim_a_sim
                k = (j-1)*dim_t_sim + [1:dim_t_sim+1:min(dim_t_sim^2,dim_t_sim*dim_a_sim - (j-1)*dim_t_sim)]';
                expI(k) = exp(-[0;cumsum(0.5*dt*(  alpha_i(k(1:end-1)) + alpha_i(k(2:end)) ...
                                                 +  beta_i(k(1:end-1)) +  beta_i(k(2:end))))]);
                if strcmp(options.gradient,'yes')
                    K = repmat(k(:),1,dim_theta);
                    L = repmat(1:dim_theta,length(k),1);
                    KL = K+(L-1)*dim_t_sim*dim_a_sim;
                    dexpIdtheta(KL) = expI(K).*(-[zeros(1,dim_theta);cumsum(0.5*dt*(  dalpha_idtheta(KL(1:end-1,:)) + dalpha_idtheta(KL(2:end,:)) ...
                                                                                    +  dbeta_idtheta(KL(1:end-1,:)) +  dbeta_idtheta(KL(2:end,:))),1)]);
                end
            end
            

        % Age-dependent division and/or death rates
        case 'age-dependent'
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{1}(1,a_sim,theta);
            beta_i  =  M.beta_fun{1}(1,a_sim,theta);
            if strcmp(options.gradient,'yes')
                dalpha_idtheta = zeros(dim_a_sim,dim_theta);
                dbeta_idtheta  = zeros(dim_a_sim,dim_theta);
                for l = 1:dim_theta
                    dalpha_idtheta(:,l) = M.dalphadtheta_fun{1,l}(1,a_sim,theta);
                    dbeta_idtheta(:,l)  = M.dbetadtheta_fun{1,l}(1,a_sim,theta);
                end
            end
            % Computation of exponential
            I = 0.5*dt*(  alpha_i(1:end-1) + alpha_i(2:end) ...
                         + beta_i(1:end-1) +  beta_i(2:end) );
            expI = zeros(dim_t_sim,dim_a_sim);
            if strcmp(options.gradient,'yes')
                dIdtheta = 0.5*dt*(  dalpha_idtheta(1:end-1,:) + dalpha_idtheta(2:end,:) ...
                                   +  dbeta_idtheta(1:end-1,:) +  dbeta_idtheta(2:end,:) );
                dexpIdtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
            end
            k = min(1:dim_a_sim,dim_t_sim);
            for j = 1:dim_a_sim
                expI(1:k(j),j) = exp(-[0;cumsum(I(k(j)-1:-1:1))]);
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        dexpIdtheta(1:k(j),j,l) = exp(-[0;cumsum(I(k(j)-1:-1:1))]).*(-[0;cumsum(dIdtheta(k(j)-1:-1:1,l))]);
                    end
                end
            end

        % Time-dependent division and/or death rates
        case 'time-dependent'
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{1}(t_sim,1,theta);
            beta_i  =  M.beta_fun{1}(t_sim,1,theta);
            if strcmp(options.gradient,'yes')
                dalpha_idtheta = zeros(dim_t_sim,dim_theta);
                dbeta_idtheta  = zeros(dim_t_sim,dim_theta);
                for l = 1:dim_theta
                    dalpha_idtheta(:,l) = M.dalphadtheta_fun{1,l}(t_sim,1,theta);
                    dbeta_idtheta(:,l)  =  M.dbetadtheta_fun{1,l}(t_sim,1,theta);
                end
            end
            % Computation of exponential
            I = [0 ; 0.5*dt*cumsum(  alpha_i(1:end-1) + alpha_i(2:end) ...
                                   +  beta_i(1:end-1) +  beta_i(2:end) ) ];
            expI = zeros(dim_t_sim,dim_a_sim);
            if strcmp(options.gradient,'yes')
                dIdtheta = [zeros(dim_theta,1) ; 0.5*dt*cumsum(  dalpha_idtheta(1:end-1,:) + dalpha_idtheta(2:end,:) ...
                                                               +  dbeta_idtheta(1:end-1,:) +  dbeta_idtheta(2:end,:) )];
                dexpIdtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
            end
            k = min(1:dim_a_sim,dim_t_sim);
            for j = 1:dim_a_sim
                expI(1:k(j),j) = exp(-I(1:k(j)));
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        dexpIdtheta(1:k(j),j,l) = exp(-I(1:k(j))).*(-dIdtheta(1:k(j),l));
                    end
                end
            end

        % Neither time- nor age-dependent division and/or death rates
        case 'constant'
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{1}(1,1,theta);
            beta_i  =  M.beta_fun{1}(1,1,theta);
            if strcmp(options.gradient,'yes')
                dalpha_idtheta = zeros(1,dim_theta);
                dbeta_idtheta  = zeros(1,dim_theta);
                for l = 1:dim_theta
                    dalpha_idtheta(l) = M.dalphadtheta_fun{1,l}(1,1,theta);
                    dbeta_idtheta(l)  =  M.dbetadtheta_fun{1,l}(1,1,theta);
                end
            end
            % Computation of exponential
            expI = triu(repmat(exp(-t_sim*(alpha_i + beta_i)),1,dim_a_sim));
            if strcmp(options.gradient,'yes')
                dexpIdtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                for l = 1:dim_theta
                    dexpIdtheta(:,:,l) = expI.*triu(repmat(-t_sim*(dalpha_idtheta(l) + dbeta_idtheta(l)),1,dim_a_sim));
                end
            end
    end

    %% Computation of age-distribution
    N_a0 = M.IC.na0.int_fun(theta);
    TA = bsxfun(@plus,-t_sim,a_sim');
    p_a0 = lognpdf(TA,M.IC.na0.mu_fun(theta),M.IC.na0.sigma_fun(theta));
    n_ai = N_a0*p_a0.*expI;
    if strcmp(options.gradient,'yes')
        dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
        dN_a0dtheta = M.IC.na0.dintdtheta_fun(theta);
        mu_a0 = M.IC.na0.mu_fun(theta);
        dmu_a0dtheta = M.IC.na0.dmudtheta_fun(theta);
        sigma_a0 = M.IC.na0.sigma_fun(theta);
        dsigmadtheta_a0 = M.IC.na0.dsigmadtheta_fun(theta);
        for l = 1:dim_theta
            dn_aidtheta(:,:,l) =   dN_a0dtheta(l)*p_a0.*expI ...
                                 + N_a0*p_a0.*(   (log(TA+1e-20) - mu_a0)   /sigma_a0^2*dmu_a0dtheta(l) ...
                                               + ((log(TA+1e-20) - mu_a0).^2/sigma_a0^3 - 1/sigma_a0)*dsigmadtheta_a0(l)).*expI ...
                                 + N_a0*p_a0.*dexpIdtheta(:,:,l);
        end
    end

    %% Visualization of age-distribution
    if strcmp(options.plot,'on');
        figure;
        surf(t_sim,a_sim,n_ai'); shading interp;
        xlabel('t');
        ylabel('a');
        title('n(a,0|t)');
    end

    %% Interpolation and integration
    % Age distribution
    if ~strcmp(type,'estimation') && (dim_a >= 1)
        Sim.n_ai(:,:,1) = Tr_t*n_ai*Tr_a;
        if strcmp(options.gradient,'yes')
            for l = 1:dim_theta
                Sim.dn_aidtheta(:,:,1,l) = Tr_t*dn_aidtheta(:,:,l)*Tr_a;
            end
        end
    end
    % Cell number
    Sim.N_i(:,1) = sum((Tr_t*n_ai).*Tr_ia0,2);
    if strcmp(options.gradient,'yes')
        for l = 1:dim_theta
            Sim.dN_idtheta(:,1,l) = sum((Tr_t*dn_aidtheta(:,:,l)).*Tr_ia0,2);
        end
    end

    %% Boundary condition
    switch options.dependency
        % Time- and age-dependent division and/or death rates
        case 'time- and age-dependent'
            b_i = 2*sum((alpha_i.*n_ai).*Tr_ia_sim0,2);
            if strcmp(options.gradient,'yes')
                for l = 1:dim_theta
                    db_idtheta(:,l) = 2*sum((  dalpha_idtheta(:,:,l).*n_ai ...
                                             + alpha_i.*dn_aidtheta(:,:,l)).*Tr_ia_sim0,2);
                end
            end
        case 'age-dependent'
            b_i = 2*sum(bsxfun(@times,n_ai,alpha_i').*Tr_ia_sim0,2);
            if strcmp(options.gradient,'yes')
                for l = 1:dim_theta
                    db_idtheta(:,l) = 2*sum((  bsxfun(@times,n_ai,dalpha_idtheta(:,l)') ...
                                             + bsxfun(@times,dn_aidtheta(:,:,l),alpha_i')).*Tr_ia_sim0,2);
                end
            end
        case 'time-dependent'
            b_i = 2*alpha_i.*sum(n_ai.*Tr_ia_sim0,2);
            if strcmp(options.gradient,'yes')
                for l = 1:dim_theta
                    db_idtheta(:,l) =   2*dalpha_idtheta(:,l).*sum(n_ai.*Tr_ia_sim0,2) ...
                                      + 2*alpha_i.*sum(dn_aidtheta(:,:,l).*Tr_ia_sim0,2);
                end
            end
        case 'constant'
            b_i = 2*alpha_i*sum(n_ai.*Tr_ia_sim0,2);
            if strcmp(options.gradient,'yes')
                for l = 1:dim_theta
                    db_idtheta(:,l) =   2*dalpha_idtheta(l)*sum(n_ai.*Tr_ia_sim0,2) ...
                                      + 2*alpha_i*sum(dn_aidtheta(:,:,l).*Tr_ia_sim0,2);
                end
            end
    end
    output = b_i;
    if strcmp(options.gradient,'yes')
        doutput = db_idtheta;
    end
    

else % b) Solution of DASP model for i = 0
    %% Computation of exponential
    % Evaluation of alpha and beta
    alpha_i = M.alpha_fun{1}(t_sim,t_sim,theta);
    beta_i  =  M.beta_fun{1}(t_sim,t_sim,theta);
    if strcmp(options.gradient,'yes')
        for l = 1:dim_theta 
            dalpha_idtheta(:,l) = M.dalphadtheta_fun{1,l}(t_sim,t_sim,theta);
            dbeta_idtheta(:,l)  = M.dbetadtheta_fun{1,l}(t_sim,t_sim,theta);
        end
    end
    % Computation of exponential
    expI = exp(-[0;cumsum(0.5*dt*(  alpha_i(1:end-1) + alpha_i(2:end) ...
                                  +  beta_i(1:end-1) +  beta_i(2:end)))]);
    if strcmp(options.gradient,'yes')
        dexpIdtheta = bsxfun(@times,expI,...
            -[zeros(1,dim_theta);cumsum(0.5*dt*(  dalpha_idtheta(1:end-1,:) + dalpha_idtheta(2:end,:) ...
                                        +  dbeta_idtheta(1:end-1,:) +  dbeta_idtheta(2:end,:)),1)]);
    end

    %% Computation of age-distribution
    k_max = min(dim_t_sim^2,dim_t_sim*dim_a_sim);
    n_ai = sparse(zeros(dim_t_sim,dim_a_sim));
    n_ai(1:dim_t_sim+1:k_max) = 2*dt*M.IC.na0.int_fun(theta)*expI;
    if strcmp(options.gradient,'yes')
        dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
        dintdtheta = M.IC.na0.dintdtheta_fun(theta);
        for l = 1:dim_theta
            dn_aidtheta(dim_t_sim*dim_a_sim*(l-1)+[1:dim_t_sim+1:k_max]) ...
                =   2*dt*dintdtheta(l)*expI ...
                  + 2*dt*M.IC.na0.int_fun(theta)*dexpIdtheta(:,l);
        end
    end
    
    %% Visualization of age-distribution
    if strcmp(options.plot,'on');
        figure;
        T = repmat(t_sim ,1,dim_a_sim);
        A = repmat(a_sim',dim_t_sim,1);
        n_ai = full(n_ai);
        ind = find(n_ai(:) > 0);
        for j = 1:length(ind)
            plot3(T(ind(j))*[1,1],A(ind(j))*[1,1],[0 n_ai(ind(j))]/(2*dt),'b-'); hold on;
            plot3(T(ind(j)),A(ind(j)),n_ai(ind(j))/(2*dt),'bo');
        end
        xlabel('t');
        ylabel('a');
        title('n(a,0|t)');
        grid on;
    end

    %% Interpolation and integration
    % Age distribution
    if ~strcmp(type,'estimation') && (dim_a >= 1)
        Sim.n_ai(:,:,1) = Tr_t*n_ai*Tr_a;
        if strcmp(options.gradient,'yes')
            for l = 1:dim_theta
                Sim.dn_aidtheta(:,:,1,l) = Tr_t*dn_aidtheta(:,:,l)*Tr_a;
            end
        end
    end
    % Cell number
    Sim.N_i(:,1) = M.IC.na0.int_fun(theta)*(Tr_t*expI);
    if strcmp(options.gradient,'yes')
        Sim.dN_idtheta(:,1,:) =   (Tr_t*expI)*dintdtheta(:)' ...
                                + M.IC.na0.int_fun(theta)*(Tr_t*dexpIdtheta);
    end

    %% Boundary condition
    b_i = 2*M.IC.na0.int_fun(theta)*alpha_i.*expI;
    if strcmp(options.gradient,'yes')
        db_idtheta =   2*M.IC.na0.int_fun(theta)*(  bsxfun(@times,dalpha_idtheta,expI) + bsxfun(@times,alpha_i,dexpIdtheta)) ...
                     + 2*bsxfun(@times,alpha_i,expI)*dintdtheta(:)';
    end

    
    output = b_i;
    if strcmp(options.gradient,'yes')
        doutput = db_idtheta;
    end

end

%% SOLUTION OF DASP MODEL FOR i >= 1
% Loop: generations
for i = 2:S 
    %% Computation of state
    % exp(-\int_0^infty \alpha_i(\tilde{a}-a+t,\tilde{a}) d\tilde{a}
    switch options.dependency
        % Time- and age-dependent division and/or death rates
        case 'time- and age-dependent'
            % Birth rate
            if i >= 3
                b_i = 2*sum((alpha_i.*n_ai).*Tr_ia_sim,2); 
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        if max(l == M.dalphadtheta_nonzero{i-1})
                            db_idtheta(:,l) = 2*sum((  bsxfun(@times,dn_aidtheta(:,:,l),alpha_i) ...
                                                     + bsxfun(@times,n_ai,dalpha_idtheta(:,:,l))).*Tr_ia_sim,2);
                        else
                            db_idtheta(:,l) = 2*sum(bsxfun(@times,dn_aidtheta(:,:,l),alpha_i).*Tr_ia_sim,2);
                        end
                    end
                end
            end
            % Definition of grid
            if i == 2
                T = repmat(t_sim ,1,dim_a_sim);
                A = repmat(a_sim',dim_t_sim,1);
            end
            % Evaluation of alpha and beta
            alpha_i = reshape(M.alpha_fun{i}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
            beta_i  = reshape( M.beta_fun{i}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
            if strcmp(options.gradient,'yes')
                if i == 2
                    dalpha_idtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                    dbeta_idtheta  = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
                for l = 1:dim_theta
                    dalpha_idtheta(:,:,l) = reshape(M.dalphadtheta_fun{i,l}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
                    dbeta_idtheta(:,:,l)  = reshape(M.dbetadtheta_fun{i,l}(T(:),A(:),theta),[dim_t_sim,dim_a_sim]);
                end
            end
            % Computation of exponential
            if i == 2
                expI = zeros(dim_t_sim,dim_a_sim);
                if strcmp(options.gradient,'yes')
                    dexpIdtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
            end
            for j = 1:dim_t_sim
                k = j + [0:dim_t_sim-j]'*(dim_t_sim+1);
                expI(k) = exp(-[ 0 ; cumsum(0.5*dt*(  alpha_i(k(1:end-1)) + alpha_i(k(2:end)) ...
                                                    +  beta_i(k(1:end-1)) +  beta_i(k(2:end)) ) ) ]);
                if strcmp(options.gradient,'yes')
                    K = repmat(k(:),1,dim_theta);
                    L = repmat(1:dim_theta,length(k),1);
                    KL = K+(L-1)*dim_t_sim*dim_a_sim;
                    dexpIdtheta(KL) = expI(K).*(-[zeros(1,dim_theta);cumsum(0.5*dt*(  dalpha_idtheta(KL(1:end-1,:)) + dalpha_idtheta(KL(2:end,:)) ...
                                                                                    +  dbeta_idtheta(KL(1:end-1,:)) +  dbeta_idtheta(KL(2:end,:))),1)]);
                end
            end
            % Computation of state
            if i == 2
                n_ai = zeros(dim_t_sim,dim_a_sim);
                if strcmp(options.gradient,'yes')
                    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
            end
            for j = 1:dim_t_sim
                k = j + [0:dim_t_sim-j]'*(dim_t_sim+1);
                n_ai(k) = b_i(j)*expI(k);
                if strcmp(options.gradient,'yes')
                    K = repmat(k(:),1,dim_theta);
                    L = repmat(1:dim_theta,length(k),1);
                    KL = K+(L-1)*dim_t_sim*dim_a_sim;
                    dn_aidtheta(KL) =   b_i(j)*dexpIdtheta(KL) + db_idtheta(j+(L-1)*dim_t_sim).*expI(K);
                end
            end 

                
        % Age-dependent division and/or death rates
        case 'age-dependent'
            % Birth rate
            if i >= 3
                b_i = 2*sum(bsxfun(@times,n_ai,alpha_i').*Tr_ia_sim,2);
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        if max(l == M.dalphadtheta_nonzero{i-1})
                            db_idtheta(:,l) = 2*sum((  bsxfun(@times,dn_aidtheta(:,:,l),alpha_i') ...
                                                     + bsxfun(@times,n_ai,dalpha_idtheta(:,l)')).*Tr_ia_sim,2);
                        else
                            db_idtheta(:,l) = 2*sum(bsxfun(@times,dn_aidtheta(:,:,l),alpha_i').*Tr_ia_sim,2);
                        end
                    end
                end
            end
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{i}(1,a_sim,theta);
            beta_i  =  M.beta_fun{i}(1,a_sim,theta);
            if strcmp(options.gradient,'yes')
                if i == 2
                    dalpha_idtheta = zeros(dim_a_sim,dim_theta);
                    dbeta_idtheta  = zeros(dim_a_sim,dim_theta);
                end
                for l = 1:dim_theta
                    dalpha_idtheta(:,l) = M.dalphadtheta_fun{i,l}(1,a_sim,theta);
                    dbeta_idtheta(:,l)  = M.dbetadtheta_fun{i,l}(1,a_sim,theta);
                end
            end
            % Computation of exponential
            expI = exp(-[ 0 ; cumsum(0.5*dt*(  alpha_i(1:end-1) + alpha_i(2:end) ...
                                              + beta_i(1:end-1) +  beta_i(2:end) ) ) ]);
            if strcmp(options.gradient,'yes')
                dexpIdtheta = bsxfun(@times,expI,...
                    -[zeros(1,dim_theta);cumsum(0.5*dt*(  dalpha_idtheta(1:end-1,:) + dalpha_idtheta(2:end,:) ...
                                                        +  dbeta_idtheta(1:end-1,:) +  dbeta_idtheta(2:end,:)),1)]);
            end
            % Computation state
            if i == 2
                n_ai = zeros(dim_t_sim,dim_a_sim);
                if strcmp(options.gradient,'yes')
                    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
            end
            k = min(1:dim_t_sim,dim_a_sim);
            for j = 1:dim_t_sim
                n_ai(j,1:k(j)) = b_i(k(j):-1:1).*expI(1:k(j));
                if strcmp(options.gradient,'yes')
                    dn_aidtheta(j,1:k(j),:) =   bsxfun(@times,db_idtheta(k(j):-1:1,:),expI(1:k(j))) ...
                                              + bsxfun(@times,b_i(k(j):-1:1),dexpIdtheta(1:k(j),:));                   
                end
            end

        % Time-dependent division and/or death rates
        case 'time-dependent'
            % Birth rate
            if i >= 3
                b_i = 2*alpha_i.*sum(n_ai.*Tr_ia_sim,2);
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        if max(l == M.dalphadtheta_nonzero{i-1})
                            db_idtheta(:,l) =   2*alpha_i            .*sum(dn_aidtheta(:,:,l).*Tr_ia_sim,2) ...
                                              + 2*dalpha_idtheta(:,l).*sum(n_ai.*Tr_ia_sim,2);
                        else
                            db_idtheta(:,l) =   2*alpha_i            .*sum(dn_aidtheta(:,:,l).*Tr_ia_sim,2);
                        end
                    end
                end
            end
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{i}(t_sim,1,theta);
            beta_i  =  M.beta_fun{i}(t_sim,1,theta);
            if strcmp(options.gradient,'yes')
                if i == 2
                    dalpha_idtheta = zeros(dim_t_sim,dim_theta);
                    dbeta_idtheta  = zeros(dim_t_sim,dim_theta);
                end
                for l = 1:dim_theta
                    dalpha_idtheta(:,l) = M.dalphadtheta_fun{i,l}(t_sim,1,theta);
                    dbeta_idtheta(:,l)  =  M.dbetadtheta_fun{i,l}(t_sim,1,theta);
                end
            end
            % Computation of exponential
            expI = exp(-[ 0 ; cumsum(0.5*dt*(  alpha_i(1:end-1) + alpha_i(2:end) ...
                                              + beta_i(1:end-1) +  beta_i(2:end) ) ) ]);
            if strcmp(options.gradient,'yes')
                dexpIdtheta = bsxfun(@times,expI,...
                    -[zeros(1,dim_theta);cumsum(0.5*dt*(  dalpha_idtheta(1:end-1,:) + dalpha_idtheta(2:end,:) ...
                                                        +  dbeta_idtheta(1:end-1,:) +  dbeta_idtheta(2:end,:)),1)]);
            end
            % Computation of state
            if i == 2
                n_ai = zeros(dim_t_sim,dim_a_sim);
                if strcmp(options.gradient,'yes')
                    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
            end
            k = min(1:dim_t_sim,dim_a_sim);
            for j = 1:dim_t_sim
                n_ai(j,1:k(j)) = b_i(k(j):-1:1).*expI(1:k(j));
                if strcmp(options.gradient,'yes')
                    dn_aidtheta(j,1:k(j),:) =   bsxfun(@times,db_idtheta(k(j):-1:1,:),expI(1:k(j))) ...
                                              + bsxfun(@times,b_i(k(j):-1:1),dexpIdtheta(1:k(j),:));                   
                end
            end
            
        
        % Neither time- nor age-dependent division and/or death rates
        case 'constant'
            % Birth rate
            if i >= 3
                b_i = 2*alpha_i*sum(n_ai.*Tr_ia_sim,2);
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        if max(l == M.dalphadtheta_nonzero{i-1})
                            db_idtheta(:,l) =   2*alpha_i          *sum(dn_aidtheta(:,:,l).*Tr_ia_sim,2) ...
                                              + 2*dalpha_idtheta(l)*sum(n_ai.*Tr_ia_sim,2);
                        else
                            db_idtheta(:,l) =   2*alpha_i          *sum(dn_aidtheta(:,:,l).*Tr_ia_sim,2);
                        end
                    end
                end
            end
            % Evaluation of alpha and beta
            alpha_i = M.alpha_fun{i}(1,1,theta);
            beta_i  =  M.beta_fun{i}(1,1,theta);
            if strcmp(options.gradient,'yes')
                if i == 2
                    dalpha_idtheta = zeros(1,dim_theta);
                    dbeta_idtheta  = zeros(1,dim_theta);
                end
                for l = 1:dim_theta
                    dalpha_idtheta(l) = M.dalphadtheta_fun{i,l}(1,1,theta);
                    dbeta_idtheta(l)  =  M.dbetadtheta_fun{i,l}(1,1,theta);
                end
            end
            % Computation of exponential
            expI = exp(-t_sim*(alpha_i + beta_i));
            if strcmp(options.gradient,'yes')
                dexpIdtheta = (-t_sim.*exp(-t_sim*(alpha_i + beta_i)))*(dalpha_idtheta + dbeta_idtheta);
            end
            % Computation state
            if i == 2
                n_ai = zeros(dim_t_sim,dim_a_sim);
                if strcmp(options.gradient,'yes')
                    dn_aidtheta = zeros(dim_t_sim,dim_a_sim,dim_theta);
                end
            end
            k = min(1:dim_t_sim,dim_a_sim);
            for j = 1:dim_t_sim
                n_ai(j,1:k(j)) = b_i(k(j):-1:1).*expI(1:k(j));
                if strcmp(options.gradient,'yes')
                    dn_aidtheta(j,1:k(j),:) =   bsxfun(@times,db_idtheta(k(j):-1:1,:),expI(1:k(j))) ...
                                              + bsxfun(@times,b_i(k(j):-1:1),dexpIdtheta(1:k(j),:));                   
                end
            end
    end

    %% Visualization of age_distribution in subpopulation
    if strcmp(options.plot,'on');
        figure;
        surf(t_sim,a_sim,n_ai'); shading interp;
        xlabel('t');
        ylabel('a');
        title(['n(a,' num2str(i-1) '|t)']);
    end

    %% Interpolation and integration
    if ~strcmp(type,'estimation')
        % Age distribution
        Sim.n_ai(:,:,i) = Tr_t*n_ai*Tr_a;
        if strcmp(options.gradient,'yes')
            for l = 1:dim_theta
                Sim.dn_aidtheta(:,:,i,l) = Tr_t*dn_aidtheta(:,:,l)*Tr_a;
            end
        end
        % Birth rate
        Sim.b_i(:,i) = Tr_t*b_i;
    end
    Sim.N_i(:,i) = sum((Tr_t*n_ai).*Tr_ia,2);
    if strcmp(options.gradient,'yes')
        for l = 1:dim_theta
            Sim.dN_idtheta(:,i,l) = sum((Tr_t*dn_aidtheta(:,:,l)).*Tr_ia,2);
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLUTION OF SET OF LSP MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The solution of the PDE for label distribution 'p(x|i,t)' is computed 
% analytically using the method of characteristics.
gamma = M.gamma_fun(theta);
k_deg = M.degradation.k_fun(theta);
c_deg = M.degradation.c_fun(theta);
mu_noise = M.noise.mu_fun(theta);
sigma_noise = M.noise.sigma_fun(theta);
px0_int = M.IC.px0.int_fun(theta);
px0_mu  = M.IC.px0.mu_fun(theta);
px0_sigma = M.IC.px0.sigma_fun(theta);
Spx0_int = sum(px0_int);
if strcmp(options.gradient,'yes')
    dgammadtheta = M.dgammadtheta_fun(theta);         dgammadtheta = dgammadtheta(:)';
    dk_degdtheta = M.degradation.dkdtheta_fun(theta); dk_degdtheta = dk_degdtheta(:)';
    dc_degdtheta = M.degradation.dcdtheta_fun(theta); dc_degdtheta = dc_degdtheta(:)';
    dmu_noisedtheta = M.noise.dmudtheta_fun(theta);   dmu_noisedtheta = dmu_noisedtheta(:)';
    dsigma_noisedtheta = M.noise.dsigmadtheta_fun(theta); dsigma_noisedtheta = dsigma_noisedtheta(:)';
    dSpx0_intdtheta = zeros(1,dim_theta);
    for k = 1:length(M.IC.px0.dintdtheta_fun)
        dpx0_intdtheta{k} = M.IC.px0.dintdtheta_fun{k}(theta); dpx0_intdtheta{k} = dpx0_intdtheta{k}(:)';
        dSpx0_intdtheta = dSpx0_intdtheta + dpx0_intdtheta{k}; 
        dpx0_mudtheta{k}  = M.IC.px0.dmudtheta_fun{k}(theta);  dpx0_mudtheta{k} = dpx0_mudtheta{k}(:)';
        dpx0_sigmadtheta{k} = M.IC.px0.dsigmadtheta_fun{k}(theta); dpx0_sigmadtheta{k} = dpx0_sigmadtheta{k}(:)';
    end
end

% Loop: generations
for i = 0:S-1
    % Loop: time points
    for j = 1:dim_t
        % Loop: components of initial distribution
        for k = 1:length(M.IC.px0.int_fun(theta))
            % Type of degradation
            if c_deg == 0
                % Exponential degradation: nu(t,x) = -k*x
                muit_k = - i*log(gamma) - k_deg*t(j) + px0_mu(k);
                sigmait_k = px0_sigma(k);
                if strcmp(options.gradient,'yes')
                    dmuit_kdtheta = - i/gamma*dgammadtheta - dk_degdtheta*t(j) + dpx0_mudtheta{k};
                    dsigmait_kdtheta = dpx0_sigmadtheta{k};
                end
            else
                % Degradation with reducing rate: nu(t,x) = -k*exp(-c*t)*x
                %  -> Gompertz decy process
                muit_k = - i*log(gamma) - k_deg/c_deg*(1-exp(-c_deg*t(j))) + px0_mu(k);
                sigmait_k = px0_sigma(k);
                if strcmp(options.gradient,'yes')
                    dmuit_kdtheta = - i/gamma*dgammadtheta ...
                        - (dk_degdtheta/c_deg - k_deg/c_deg^2*dc_degdtheta)*(1-exp(-c_deg*t(j))) ...
                        - k_deg/c_deg*exp(-c_deg*t(j))*t(j)*dc_degdtheta  ...
                        + dpx0_mudtheta{k};
                    dsigmait_kdtheta = dpx0_sigmadtheta{k};
                end
            end
            % Summation of contribution of subpopulations
            if ~strcmp(type,'estimation')
                lognpdfx = lognpdf(x',muit_k,sigmait_k);
                Sim.p_xi(j,:,i+1) = Sim.p_xi(j,:,i+1) ...
                    + px0_int(k)/Spx0_int*lognpdfx;
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        Sim.dp_xidtheta(j,:,i+1,l) = Sim.dp_xidtheta(j,:,i+1,l) ...
                            + (dpx0_intdtheta{k}(l)/Spx0_int-px0_int(k)/Spx0_int*dSpx0_intdtheta(l))*lognpdfx ...
                            + px0_int(k)/Spx0_int*lognpdfx.*(...
                                        ((log(x')-muit_k)/sigmait_k^2)*dmuit_kdtheta(l)   ...
                                      + ((log(x')-muit_k).^2/sigmait_k^3 - 1/sigmait_k)*dsigmait_kdtheta(l));
                    end
                end
            end
            
            % Computation of noise corrupted measured label distribution
            % in subpopulations
            if strcmp(options.noise.flag,'yes')
                % Mean and variance for label distribution
                mean1 = exp(muit_k)*exp(sigmait_k^2/2);
                var1 = exp(2*muit_k)*exp(sigmait_k^2)*(exp(sigmait_k^2)-1);
                if strcmp(options.gradient,'yes')
                    dmean1dtheta = exp(muit_k)*exp(sigmait_k^2/2) * (dmuit_kdtheta + sigmait_k*dsigmait_kdtheta);
                    dvar1dtheta =   exp(2*muit_k)*exp(sigmait_k^2)*(2*exp(sigmait_k^2)*sigmait_k*dsigmait_kdtheta) ...
                                  + exp(2*muit_k)*(2*exp(sigmait_k^2)*sigmait_k*dsigmait_kdtheta)*(exp(sigmait_k^2)-1) ...
                                  + (2*exp(2*muit_k)*dmuit_kdtheta)*exp(sigmait_k^2)*(exp(sigmait_k^2)-1);
                end
                % Mean and variance for measurement noise
                mean2 = exp(mu_noise)*exp(sigma_noise^2/2);
                var2 = exp(2*mu_noise)*exp(sigma_noise^2)*(exp(sigma_noise^2)-1);
                if strcmp(options.gradient,'yes')
                    dmean2dtheta = exp(mu_noise)*exp(sigma_noise^2/2)*(dmu_noisedtheta + sigma_noise*dsigma_noisedtheta);
                    dvar2dtheta =   exp(2*mu_noise)*exp(sigma_noise^2)*(2*sigma_noise*dsigma_noisedtheta) ...
                                  + exp(2*mu_noise)*(2*exp(sigma_noise^2)*sigma_noise*dsigma_noisedtheta)*(exp(sigma_noise^2)-1) ...
                                  + (2*exp(2*mu_noise)*dmu_noisedtheta)*exp(sigma_noise^2)*(exp(sigma_noise^2)-1);
                end
                % Mean and variance of convoluted distribution
                mean = mean1 + mean2;
                var = var1 + var2;
                if strcmp(options.gradient,'yes')
                    dmeandtheta = dmean1dtheta + dmean2dtheta;
                    dvardtheta = dvar1dtheta + dvar2dtheta;
                end
                % mu and sigma of convoluted distribution
                sigma = sqrt(log(var/mean^2 + 1));
                mu = log(mean) - sigma^2/2;
                if strcmp(options.gradient,'yes')
                    dsigmadtheta = 0.5*(log(var/mean^2 + 1))^(-1/2) / (var/mean^2 + 1) * (dvardtheta/mean^2 - 2*var*dmeandtheta/mean^3);
                    dmudtheta = dmeandtheta/mean - sigma*dsigmadtheta;
                end
                % Summation of contribution of subpopulations
                lognpdfx = lognpdf(x',mu,sigma);
                Sim.p_yi(j,:,i+1) = Sim.p_yi(j,:,i+1) + ...
                    px0_int(k)*lognpdfx;
                if strcmp(options.gradient,'yes')
                    for l = 1:dim_theta
                        Sim.dp_yidtheta(j,:,i+1,l) = Sim.dp_yidtheta(j,:,i+1,l) ...
                            + (dpx0_intdtheta{k}(l)/Spx0_int-px0_int(k)/Spx0_int*dSpx0_intdtheta(l))*lognpdfx ...
                            + px0_int(k)/Spx0_int*lognpdfx.*(...
                                        ((log(x')-mu)/sigma^2)*dmudtheta(l)   ...
                                      + ((log(x')-mu).^2/sigma^3 - 1/sigma)*dsigmadtheta(l));
                    end
                end
            end
        end
    end
end
% Noise-free case:
if ~strcmp(options.noise.flag,'yes')
    Sim.p_yi = Sim.p_xi;
    Sim.dp_yidtheta = Sim.dp_xidtheta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATION OF OVERALL RESPONSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cell number
Sim.N = sum(Sim.N_i,2);
if strcmp(options.gradient,'yes')
    Sim.dNdtheta = squeeze(sum(Sim.dN_idtheta,2));
end
    
%% Label and fluorescence distribution
if ~strcmp(type,'estimation')
    % Subpopulation
    Sim.n_xi = permute(repmat(Sim.N_i,[1,1,dim_x]),[1,3,2]).*Sim.p_xi;
    Sim.n_xi = squeeze(Sim.n_xi);
    Sim.n_yi = permute(repmat(Sim.N_i,[1,1,dim_x]),[1,3,2]).*Sim.p_yi;
    Sim.n_yi = squeeze(Sim.n_yi);
    % Overall population
    Sim.n_x = sum(Sim.n_xi,3);
    Sim.n_y = sum(Sim.n_yi,3);
    % Normalized overall label distribution
    Sim.p_y = bsxfun(@times,Sim.n_y,1./Sim.N);
else
    % Overall fluorescence distribution
    Sim.n_y = sum(permute(repmat(Sim.N_i,[1,1,dim_x]),[1,3,2]).*Sim.p_yi,3);
    Sim.p_y = sum(permute(repmat(bsxfun(@rdivide,Sim.N_i,Sim.N),[1,1,dim_x]),[1,3,2]).*Sim.p_yi,3);
end
% Gradients
if strcmp(options.gradient,'yes')
    for l = 1:dim_theta
        Sim.dn_ydtheta(:,:,l) =   sum(permute(repmat(Sim.dN_idtheta(:,:,l),[1,1,dim_x]),[1,3,2]).*Sim.p_yi,3) ...
                                + sum(permute(repmat(Sim.N_i,[1,1,dim_x]),[1,3,2]).*Sim.dp_yidtheta(:,:,:,l),3);
        Sim.dp_ydtheta(:,:,l) =   sum(permute(repmat(bsxfun(@rdivide,Sim.dN_idtheta(:,:,l),Sim.N),[1,1,dim_x]),[1,3,2]).*Sim.p_yi,3) ...
                                - sum(permute(repmat(bsxfun(@times,bsxfun(@rdivide,Sim.N_i,(Sim.N.^2)),Sim.dNdtheta(:,l)),[1,1,dim_x]),[1,3,2]).*Sim.p_yi,3) ...
                                + sum(permute(repmat(bsxfun(@rdivide,Sim.N_i,Sim.N),[1,1,dim_x]),[1,3,2]).*Sim.dp_yidtheta(:,:,:,l),3);
    end
end

%% Average division number
if ~strcmp(type,'estimation')
    Sim.ave_i = sum(Sim.N_i*diag([0:S-1]),2)./Sim.N;
end

%%
if ~exist('output')
    output = [];
    doutput = [];
end