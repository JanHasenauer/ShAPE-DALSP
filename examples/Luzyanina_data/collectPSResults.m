n_max = 250;

%% Setting
switch M_number
    case 1
        model_1__alpha__beta__a0_delta;
    case 2
        model_2__alpha_a__beta__a0_delta;
    case 16
        model_16__alpha_ia__beta_ia__a0_delta;
        
end
Options_PS.Size    = 1e1*parameters.number;
Options_PS.MaxIter = 1e3*parameters.number;
Options_PS.MaxObj  = 1e3*parameters.number;

%% Collect data
parameters.MS.par = nan(parameters.number,n_max);
parameters.MS.logPost = nan(n_max,1);
parameters.MS.n_objfun = nan(n_max,1);
parameters.MS.n_iter = nan(n_max,1);
parameters.MS.t_cpu = nan(n_max,1);

for i = 1:n_max
    try 
        path = [ './project/results/results'  num2str(M_number,'%d') '__PS__Size_' num2str(Options_PS.Size,'%d') '__MaxObj_' num2str(Options_PS.MaxObj,'%d') '__Start_' num2str(i)];
        parameters.MS.par(:,i) = csvread([path '_theta']);
        parameters.MS.logPost(i) = -csvread([path '_J']);
        parameters.MS.n_objfun(i) = csvread([path '_ObjFunCounter']);
        parameters.MS.n_iter(i) = csvread([path '_IterCounter']);
        parameters.MS.t_cpu(i) = csvread([path '_tcpu']);
    end
end
parameters = sortMultiStarts(parameters);
