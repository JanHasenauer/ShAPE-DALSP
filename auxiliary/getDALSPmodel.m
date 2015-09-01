% getDALSPmodel converts a symbolically defined DALSP models in a series 
%   of numerical functions. The evaluation of these functions provides 
%   for examples the rate of cell division and its derivative with respect 
%   to the model parameters.

function M = getDALSPmodel(M,parameters)

%% INITIALIZATION
tdep = 0;
adep = 0;

%% DIVISION AND DEATH RATES
% FUNCTIONAL EXPRESSION FOR MODEL QUANTITIES
% alpha
for i = 1:M.S
    if M.alpha{i} == 0
        str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
    else
        % Differentiation
        Dt = simplify(diff(M.alpha{i},M.time));
        Da = simplify(diff(M.alpha{i},M.age));
        % Conversion of functional experssion to string
        % constant
        if (Dt == 0) && (Da == 0)
            str = ['@(t,a,theta) ones(max(size(a),size(t)))*' char(M.alpha{i})];
        end
        % time-dependent
        if (Dt ~= 0) && (Da == 0)
            str = ['@(t,a,theta) ones(size(a))*' char(M.alpha{i})];
            tdep = 1;
        end
        % age-dependent
        if (Dt == 0) && (Da ~= 0)
            str = ['@(t,a,theta) ones(size(t))*' char(M.alpha{i})];
            adep = 1;
        end
        % time- and age-dependent
        if (Dt ~= 0) && (Da ~= 0)
            str = ['@(t,a,theta) ' char(M.alpha{i})];
            tdep = 1; adep = 1;
        end
        % Replacement of parameters
        for k = 1:parameters.number
            str = strrep(str,char(parameters.sym(k)),['theta(' num2str(k,'%d') ')']);
        end
        % Vectorization
        str = strrep(str,'*','.*');
        str = strrep(str,'/','./');
        str = strrep(str,'^','.^');
    end
    % Conversion to function handle
    M.alpha_fun{i} = eval(str);
end

% beta
for i = 1:M.S
    if M.beta{i} == 0
        str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
    else
        % Differentiation
        Dt = simplify(diff(M.beta{i},M.time));
        Da = simplify(diff(M.beta{i},M.age));
        % Conversion of functional experssion to string
        % constant
        if (Dt == 0) && (Da == 0)
            str = ['@(t,a,theta) ones(max(size(a),size(t)))*' char(M.beta{i})];
        end
        % time-dependent
        if (Dt ~= 0) && (Da == 0)
            str = ['@(t,a,theta) ones(size(a))*' char(M.beta{i})];
            tdep = 1;
        end
        % age-dependent
        if (Dt == 0) && (Da ~= 0)
            str = ['@(t,a,theta) ones(size(t))*' char(M.beta{i})];
            adep = 1;
        end
        % time- and age-dependent
        if (Dt ~= 0) && (Da ~= 0)
            str = ['@(t,a,theta) ' char(M.beta{i})];
            tdep = 1; adep = 1;
        end
        % Replacement of parameters
        for k = 1:parameters.number
            str = strrep(str,char(parameters.sym(k)),['theta(' num2str(k,'%d') ')']);
        end
        % Vectorization
        str = strrep(str,'*','.*');
        str = strrep(str,'/','./');
        str = strrep(str,'^','.^');
    end
    % Conversion to function handle
    M.beta_fun{i} = eval(str);
end

% Asignment/Check of model type
if isfield(M,'type')
    if ~strcmp(M.type,'constant')       && ~strcmp(M.type,'age-dependent') && ...
       ~strcmp(M.type,'time-dependent') && ~strcmp(M.type,'time- and age-dependent')
        error('Unknown model type (M.type). Only ''constant'', ''time-dependent'', ''age-dependent'', and ''time- and age-dependent'' are allowed.')
    end
end
% constant
if (tdep == 0) && (adep == 0)
    type = 'constant';
    if ~isfield(M,'type')
        M.type = type;
    end
end
% time-dependent
if (tdep ~= 0) && (adep == 0)
    type = 'time-dependent';
    if isfield(M,'type')
        if strcmp(M.type,'constant') || strcmp(M.type,'age-dependent')
            error('Model type (M.type) is inconsistent. The prolifertaion rates are time-dependent.')
        end
    else
        M.type = type;
    end
end
% age-dependent
if (tdep == 0) && (adep ~= 0)
    type = 'age-dependent';
    if isfield(M,'type')
        if strcmp(M.type,'constant') || strcmp(M.type,'time-dependent')
            error('Model type (M.type) is inconsistent. The prolifertaion rates are age-dependent.')
        end
    else
        M.type = type;
    end
end
% time- and age-dependent
if (tdep ~= 0) && (adep ~= 0)
    type = 'time- and age-dependent';
    if isfield(M,'type')
        if strcmp(M.type,'constant') || strcmp(M.type,'time-dependent') || strcmp(M.type,'age-dependent')
            error('Model type (M.type) is inconsistent. The prolifertaion rates are time- and age-dependent.')
        end
    else
        M.type = type;
    end
end

% FUNCTIONAL EXPRESSION FOR DERIVATIVES OF MODEL QUANTITIES
% alpha
for i = 1:M.S
    % Index set of non-zero derivatives
    M.dalphadtheta_nonzero{i} = [];
    % Loop: parameters
    for j = 1:parameters.number
        if M.alpha{i} == 0
            str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
        else
            % Differentiation with respect to parameters
            D = simplify(diff(M.alpha{i},parameters.sym(j)));
            % Check whether derivative is zero
            if D == 0
                str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
            else
                % Differentiation
                Dt = simplify(diff(D,M.time));
                Da = simplify(diff(D,M.age));
                % Conversion of functional experssion to string
                % constant
                if (Dt == 0) && (Da == 0)
                    str = ['@(t,a,theta) ones(max(size(a),size(t)))*' char(D)];
                end
                % time-dependent
                if (Dt ~= 0) && (Da == 0)
                    str = ['@(t,a,theta) ones(size(a))*' char(D)];
                end
                % age-dependent
                if (Dt == 0) && (Da ~= 0)
                    str = ['@(t,a,theta) ones(size(t))*' char(D)];
                end
                % time- and age-dependent
                if (Dt ~= 0) && (Da ~= 0)
                    str = ['@(t,a,theta) ' char(D)];
                end
                % Replacement of parameters
                for k = 1:parameters.number
                    str = strrep(str,char(parameters.sym(k)),['theta(' num2str(k,'%d') ')']);
                end
                % Vectorization
                str = strrep(str,'*','.*');
                str = strrep(str,'/','./');
                str = strrep(str,'^','.^');
                % Compensation of numerical problem with a*log(a)
                str = strrep(str,'log(a)','log(a+1e-20)');
                str = strrep(str,'log(t)','log(t+1e-20)');
                % Update index set
                M.dalphadtheta_nonzero{i} = [M.dalphadtheta_nonzero{i},j];
            end
        end
        % Conversion to function handle
        M.dalphadtheta_fun{i,j} = eval(str);
    end
end

% beta
for i = 1:M.S
    % Index set of non-zero derivatives
    M.dbetadtheta_nonzero{i} = [];
    % Loop: parameters
    for j = 1:parameters.number
        if M.beta{i} == 0
            str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
        else
            % Differentiation with respect to parameters
            D = simplify(diff(M.beta{i},parameters.sym(j)));
            % Check whether derivative is zero
            if D == 0
                str = ['@(t,a,theta) zeros(max(size(a),size(t)))'];
            else
                % Differentiation
                Dt = simplify(diff(D,M.time));
                Da = simplify(diff(D,M.age));
                % Conversion of functional experssion to string
                % constant
                if (Dt == 0) && (Da == 0)
                    str = ['@(t,a,theta) ones(max(size(a),size(t)))*' char(D)];
                end
                % time-dependent
                if (Dt ~= 0) && (Da == 0)
                    str = ['@(t,a,theta) ones(size(a))*' char(D)];
                end
                % age-dependent
                if (Dt == 0) && (Da ~= 0)
                    str = ['@(t,a,theta) ones(size(t))*' char(D)];
                end
                % time- and age-dependent
                if (Dt ~= 0) && (Da ~= 0)
                    str = ['@(t,a,theta) ' char(D)];
                end
                % Replacement of parameters
                for k = 1:parameters.number
                    str = strrep(str,char(parameters.sym(k)),['theta(' num2str(k,'%d') ')']);
                end
                % Vectorization
                str = strrep(str,'*','.*');
                str = strrep(str,'/','./');
                str = strrep(str,'^','.^');
                % Compensation of numerical problem with a*log(a)
                str = strrep(str,'log(a)','log(a+1e-20)');
                str = strrep(str,'log(t)','log(t+1e-20)');
                % Update index set
                M.dbetadtheta_nonzero{i} = [M.dbetadtheta_nonzero{i},j];
            end
        end
        % Conversion to function handle
        M.dbetadtheta_fun{i,j} = eval(str);
    end
end

%% OTHER PARAMETERS
% FUNCTIONAL EXPRESSION FOR MODEL QUANTITIES

% gamma
M.gamma_fun = sym2fun(M.gamma,parameters.sym);

% degradation
M.degradation.k_fun = sym2fun(M.degradation.k,parameters.sym);
M.degradation.c_fun = sym2fun(M.degradation.c,parameters.sym);

% noise
M.noise.mu_fun    = sym2fun(M.noise.mu   ,parameters.sym);
M.noise.sigma_fun = sym2fun(M.noise.sigma,parameters.sym);

% noise N
M.noise_N.sigma_fun = sym2fun(M.noise_N.sigma,parameters.sym);

% initial condition
M.IC.na0.int_fun   = sym2fun(M.IC.na0.int  ,parameters.sym);
try
    M.IC.na0.mu_fun    = sym2fun(M.IC.na0.mu   ,parameters.sym);    
    M.IC.na0.sigma_fun = sym2fun(M.IC.na0.sigma,parameters.sym);
catch
    M.IC.na0.type = 'delta';
end
M.IC.px0.int_fun   = sym2fun(M.IC.px0.int  ,parameters.sym);
M.IC.px0.mu_fun    = sym2fun(M.IC.px0.mu   ,parameters.sym);
M.IC.px0.sigma_fun = sym2fun(M.IC.px0.sigma,parameters.sym);


% FUNCTIONAL EXPRESSION FOR DERIVATIVES OF MODEL QUANTITIES

% gamma
M.dgammadtheta_fun = sym2fun(jacobian(M.gamma,parameters.sym),parameters.sym);

% degradation
M.degradation.dkdtheta_fun = sym2fun(jacobian(M.degradation.k,parameters.sym),parameters.sym);
M.degradation.dcdtheta_fun = sym2fun(jacobian(M.degradation.c,parameters.sym),parameters.sym);

% noise
M.noise.dmudtheta_fun    = sym2fun(jacobian(M.noise.mu   ,parameters.sym),parameters.sym);
M.noise.dsigmadtheta_fun = sym2fun(jacobian(M.noise.sigma,parameters.sym),parameters.sym);

% noise N
M.noise_N.dsigmadtheta_fun = sym2fun(jacobian(M.noise_N.sigma,parameters.sym),parameters.sym);

% initial condition
M.IC.na0.dintdtheta_fun = sym2fun(jacobian(M.IC.na0.int  ,parameters.sym),parameters.sym);
try
    M.IC.na0.dmudtheta_fun    = sym2fun(jacobian(M.IC.na0.mu   ,parameters.sym),parameters.sym);    
    M.IC.na0.dsigmadtheta_fun = sym2fun(jacobian(M.IC.na0.sigma,parameters.sym),parameters.sym);
end
for i = 1:length(M.IC.px0.int)
    M.IC.px0.dintdtheta_fun{i}   = sym2fun(jacobian(M.IC.px0.int(i)  ,parameters.sym),parameters.sym);
    M.IC.px0.dmudtheta_fun{i}    = sym2fun(jacobian(M.IC.px0.mu(i)   ,parameters.sym),parameters.sym);
    M.IC.px0.dsigmadtheta_fun{i} = sym2fun(jacobian(M.IC.px0.sigma(i),parameters.sym),parameters.sym);
end
       
