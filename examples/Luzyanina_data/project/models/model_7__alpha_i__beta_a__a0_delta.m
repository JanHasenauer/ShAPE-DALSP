%% SPECIFY MODEL PARAMETERS
% Symbolic quantities
syms a t;
% Symbolic parameters
syms log10_k_alpha_0 ...
     log10_k_alpha_1 ...
     log10_k_alpha_2 ...
     log10_k_alpha_3 ...
     log10_k_alpha_4 ...
     log10_k_alpha_5 ...
     log10_k_alpha_6 ...
     log10_K_beta    log10_n_beta ...
     log10_k_beta    ...
     log10_k_deg     log10_c_deg       ...
     mu_noise        log10_sigma_noise ...
     log10_N_0       r_x1              ...
     mu_x1           log10_sigma_x1    ...
     mu_x2           log10_sigma_x2;
% Symbolic parameter vector
parameters.sym  = [ log10_k_alpha_0 ; ...
                    log10_k_alpha_1 ; ...
                    log10_k_alpha_2 ; ...
                    log10_k_alpha_3 ; ...
                    log10_k_alpha_4 ; ...
                    log10_k_alpha_5 ; ...
                    log10_k_alpha_6 ; ...
                    log10_K_beta    ; log10_n_beta ; ...
                    log10_k_beta    ; ...
                    log10_k_deg     ; log10_c_deg        ; ...
                    mu_noise        ; log10_sigma_noise  ; ...
                    log10_N_0       ; r_x1               ; ...
                    mu_x1           ; log10_sigma_x1     ; ...
                    mu_x2           ; log10_sigma_x2 ];
% Names of parameters as strings
parameters.name = {'log_{10}(k_{\alpha,0})'; ...
                   'log_{10}(k_{\alpha,1})'; ...
                   'log_{10}(k_{\alpha,2})'; ...
                   'log_{10}(k_{\alpha,3})'; ...
                   'log_{10}(k_{\alpha,4})'; ...
                   'log_{10}(k_{\alpha,5})'; ...
                   'log_{10}(k_{\alpha,6})'; ...
                   'log_{10}(K_{\beta})'   ;'log_{10}(n_{\beta})' ; ...
                   'log_{10}(k_{\beta})'   ; ...
                   'log_{10}(k_{deg})'     ;'log_{10}(c_{deg})'       ; ...
                   '\mu_{noise}'           ;'log_{10}(\sigma_{noise})'; ...
                   'log_{10}(N_0)'         ; 'r_{x,1}'                ; ...
                   '\mu_{x,1}'             ; 'log_{10}(\sigma_{x,1})' ; ...
                   '\mu_{x,2}'             ; 'log_{10}(\sigma_{x,2})'};
% Number of parameter values
parameters.number = length(parameters.sym);
% Initial guess of parameter value
parameters.guess = [-0.1224; ...
                     0.8093; ...
                     0.4509; ...
                     0.3467; ...
                     2.0000; ...
                     2.0000; ...
                     2.0000; ...
                    -0.1603; 0.7610; ...
                    -3.0401; ...
                    -0.7293; -0.9206; ...
                     0.9396; -0.3986; ...
                     3.8108;  0.5717; ...
                     6.3758; -0.8326; ...
                     6.7665; -0.3358];
% Minimum and maximum for parameter values
parameters.min  = [  -6; ...
                     -6; ...
                     -6; ...
                     -6; ...
                     -6; ...
                     -6; ...
                     -6; ...
                     -6;  -6; ...
                     -6; ...
                     -4; -4; ...
                      0; -1; ...
                    3.5;  0; ...
                      6; -1; ...
                      6; -1];

parameters.max  = [   3; ...
                      3; ...
                      3; ...
                      3; ...
                      3; ...
                      3; ...
                      3; ...
                      2;   2; ...
                      3; ...
                      0;  0; ...
                      4;  1; ... 
                      4;  1; ...
                      8;  0; ...
                      8;  0];

%% MODEL
M.type = 'age-dependent';
M.time = t;
M.age = a;
M.name = 'model_7__alpha_i__beta_a__a0_delta';
% Number of subpopulations
M.S = 8;
% Initial condition
M.IC.na0.int = 10^log10_N_0;
M.IC.na0.type = 'delta';
M.IC.px0.int = [ r_x1;1-r_x1];
M.IC.px0.mu  = [mu_x1; mu_x2];
M.IC.px0.sigma = [10^log10_sigma_x1;10^log10_sigma_x2];
% Division rates
M.alpha{1} = 10^log10_k_alpha_0;
M.alpha{2} = 10^log10_k_alpha_1;
M.alpha{3} = 10^log10_k_alpha_2;
M.alpha{4} = 10^log10_k_alpha_3;
M.alpha{5} = 10^log10_k_alpha_4;
M.alpha{6} = 10^log10_k_alpha_5;
M.alpha{7} = 10^log10_k_alpha_6;
M.alpha{8} = 0;
% Death rates
M.beta{1} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{2} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{3} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{4} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{5} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{6} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{7} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
M.beta{8} = (10^log10_k_beta) * (a.^(10^log10_n_beta))./((10^log10_K_beta)^(10^log10_n_beta) + a.^(10^log10_n_beta));
% Remaining parameters
M.gamma  = 2;
M.degradation.k = 10^log10_k_deg;
M.degradation.c = 10^log10_c_deg;
M.noise.mu      = mu_noise;
M.noise.sigma   = 10^log10_sigma_noise;
M.noise_N.sigma = 0.01;

% Generate model
M = getDALSPmodel(M,parameters);