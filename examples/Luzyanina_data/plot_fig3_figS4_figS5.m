clear all;
close all;
clc;

%% LOAD DATA
% % main text
% load('project/results/model_1__alpha__beta__a0_delta__grad_on');        Mo{1} = M; par{1} = parameters; II(1) = 1;
% load('project/results/model_2__alpha_a__beta__a0_delta__grad_on');      Mo{2} = M; par{2} = parameters; II(2) = 2;
% load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on'); Mo{3} = M; par{3} = parameters; II(3) = 16;
% filename = 'fig3';

% % supplement
% load('project/results/model_1__alpha__beta__a0_delta__grad_on');      Mo{1} = M; par{1} = parameters; II(1) = 1;
% load('project/results/model_2__alpha_a__beta__a0_delta__grad_on');    Mo{2} = M; par{2} = parameters; II(2) = 2;
% load('project/results/model_3__alpha_i__beta__a0_delta__grad_on');    Mo{3} = M; par{3} = parameters; II(3) = 3;
% load('project/results/model_4__alpha_ia__beta__a0_delta__grad_on');   Mo{4} = M; par{4} = parameters; II(4) = 4;
% load('project/results/model_5__alpha__beta_a__a0_delta__grad_on');    Mo{5} = M; par{5} = parameters; II(5) = 5;
% load('project/results/model_6__alpha_a__beta_a__a0_delta__grad_on');  Mo{6} = M; par{6} = parameters; II(6) = 6;
% load('project/results/model_7__alpha_i__beta_a__a0_delta__grad_on');  Mo{7} = M; par{7} = parameters; II(7) = 7;
% load('project/results/model_8__alpha_ia__beta_a__a0_delta__grad_on'); Mo{8} = M; par{8} = parameters; II(8) = 8;
% filename = 'figS4';

% supplement
load('project/results/model_9__alpha__beta_i__a0_delta__grad_on');      Mo{1} = M; par{1} = parameters; II(1) = 9;
load('project/results/model_10__alpha_a__beta_i__a0_delta__grad_on');   Mo{2} = M; par{2} = parameters; II(2) = 10;
load('project/results/model_11__alpha_i__beta_i__a0_delta__grad_on');   Mo{3} = M; par{3} = parameters; II(3) = 11;
load('project/results/model_12__alpha_ia__beta_i__a0_delta__grad_on');  Mo{4} = M; par{4} = parameters; II(4) = 12;
load('project/results/model_13__alpha__beta_ia__a0_delta__grad_on');    Mo{5} = M; par{5} = parameters; II(5) = 13;
load('project/results/model_14__alpha_a__beta_ia__a0_delta__grad_on');  Mo{6} = M; par{6} = parameters; II(6) = 14;
load('project/results/model_15__alpha_i__beta_ia__a0_delta__grad_on');  Mo{7} = M; par{7} = parameters; II(7) = 15;
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on'); Mo{8} = M; par{8} = parameters; II(8) = 16;
filename = 'figS5';

%%
for i = 1:length(II)
    S{i} = CPsimulateDALSP(Mo{i},par{i}.MS.par(:,1),D.t,a_sim,...
                            [D(1).bins(1,1); D(1).bins(:,2)]',...
                            options_logL.simulation);
end

%% INITIALIZATION
%col_data = 0.7*[1,1,1];
col_data = [0.3,0.5,1];
col_sim = 0*[1,1,1];
col_sim_bounds = 0*[1,1,1];
col_bg = 0.5*[1,1,1];

x = [D.bins(:,1);D.bins(end,2)];

kmax = size(S{1}.n_y,1);

for i = 1:length(S)
    h_bg{i} = pdf('logn',x,Mo{i}.noise.mu_fun(par{i}.MS.par(:,1)),Mo{i}.noise.sigma_fun(par{i}.MS.par(:,1)));
    H_bg{i} = (h_bg{i}(1:end-1) + h_bg{i}(2:end))/2.*diff(x);

    S{i}.H = bsxfun(@times,S{i}.n_y(:,1:end-1) + S{i}.n_y(:,2:end), 0.5*diff(x(:)'));
    S{i}.h = bsxfun(@times,S{i}.p_y(:,1:end-1) + S{i}.p_y(:,2:end), 0.5*diff(x(:)'));

    lb{i} = zeros(size(S{i}.H));
    ub{i} = zeros(size(S{i}.H));
    for j = 1:size(lb{i},1)
        for k = 1:size(lb{i},2)
            b = binoinv([0.05 0.95],sum(D.cellcount(j,:)),S{i}.h(j,k));
            lb{i}(j,k) = b(1);
            ub{i}(j,k) = b(2);
        end
    end
    
    t_sim = linspace(0,D.t(end),30);
    [S2{i}] = CPsimulateDALSP(Mo{i},par{i}.MS.par(:,1),t_sim,t_sim,...
                        [D(1).bins(1,1); D(1).bins(:,2)]',...
                        options_logL.simulation);
end


%% CFSE DISTRIBUTION
figure;
ll = [];

getIsp = @(i,j) [(9*(i-1)+2+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+3+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+4+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+5+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+6+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+7+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+8+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)],...
                 (9*(i-1)+9+6)*(5*(length(D.t)+1)+2)+[(5*(j-1)+1)+(j>5):5*j+(j>5)]];

% Plot mueasurement and simulation / construct legend label
for i = 0:length(S)
for j = 1:length(D.t)
    subplot(8+9*(length(S)+1),5*(length(D.t)+1)+2,getIsp(i+1,j));
    % Background
    ind = round(1:0.5:length(x)-0.5);
    ind2 = round(0.5:0.5:length(x)-1);

    if j == 1 && i >= 1
        ll(3) = fill(x([ind,ind(end:-1:1)]),15000*[H_bg{i}(ind2);zeros(size(lb{i}(j,ind2(end:-1:1))'))],col_bg); hold on;
        stairs(x,15000*H_bg{i}([1:end,end]),'-','Color',col_bg,'linewidth',1); hold on;
    end
   
    % Measurement
    ll(1) = stairs(x,D.cellcount(j,[1:end,end]),'-','Color',col_data,'linewidth',1); hold on;

    % Simulation
    if i >= 1
    %fill(x([ind,ind(end:-1:1)]),[ub(j,ind2)';lb(j,ind2(end:-1:1))'],0.8*[1,1,1]); hold on;
    stairs(x,lb{i}(j,[1:end,end]),'-','Color',col_sim_bounds,'linewidth',1); hold on;
    stairs(x,ub{i}(j,[1:end,end]),'-','Color',col_sim_bounds,'linewidth',1); hold on;

    % Simulation
    ll(2) = stairs(x,S{i}.h(j,[1:end,end])*sum(D.cellcount(j,:)),'-','Color',col_sim,'linewidth',2); hold on;
    end
    
    % Set legend / label
    set(gca,'Xscale','log');
    %xlim(x([1,end]));
    xlim([1,x(end)]);
    ylim([0,150]);
    set(gca,'ytick',[0:50:150]);
    set(gca,'xtick',[2:9,10:10:90,100:100:900,1000],'xticklabel',[]);%{'','','','','','','10','','','','','','','','','100','','','','','','','','','1000','',''});
    if i == length(S)
        text(1,-12,'10^0','fontsize',12);
        text(10,-12,'10^1','fontsize',12);
        text(100,-12,'10^2','fontsize',12);
        text(1000,-12,'10^3','fontsize',12);
        text(4,-30,'fluorescence [UI]','fontsize',12);
    elseif i == 0
        title(['day ' num2str(j,'%d')]);
        if j == 3
            title({'{\bf CFSE distribution}',['day ' num2str(j,'%d')]});
        end
    end
    if j >= 2
        set(gca,'yticklabel',[]);
    else
        ylabel('cell count per bin [-]');
        f = 2;
        switch i
            case 0
                text(f*0.0006,146,'{\bf A}','fontsize',16);
                text(f*0.0020,145,'{\bf $(\mathcal{D})$}','Interpreter','latex','fontsize',16);
            case 1
                text(f*0.0006,146,'{\bf B}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 2
                text(f*0.0006,146,'{\bf C}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 3
                text(f*0.0006,146,'{\bf D}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 4
                text(f*0.0006,146,'{\bf E}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 5
                text(f*0.0006,146,'{\bf F}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 6
                text(f*0.0006,146,'{\bf G}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 7
                text(f*0.0006,146,'{\bf H}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
            case 8
                text(f*0.0006,146,'{\bf I}','fontsize',16);
                text(f*0.0020,145,['{\bf $(\mathcal{M}_{' num2str(II(i),'%d') '})$}'],'Interpreter','latex','fontsize',16);
        end
        if i == length(S)
            legend(ll,{'data','model','background'},'location','NorthWest');
        end
    end
end

% Cell number
j = 6;
subplot(8+9*(length(S)+1),5*(length(D.t)+1)+2,getIsp(i+1,j));
ll2(1) = plot(D.t_plot,sum(D.cellcount,2),'-','linewidth',1,'Color',col_data); hold on;
ll2(1) = plot(D.t_plot,sum(D.cellcount,2),'.','MarkerSize',20,'Color',col_data); hold on;
if i >= 1
    ll2(2) = plot(t_sim+D.t_plot(1),S2{i}.N,'-','linewidth',2,'Color',col_sim); hold on;
end

% Set legende label
xlim(D.t_plot([1,end]));
ylabel('overall cell count [-]');
set(gca,'xtick',[1:5]);
if i == length(S)
    xlabel('time [days]');
else
    set(gca,'xticklabel',[]);
end
if i == 0
    title({'{\bf Population size}',''});
end
if i == length(S)
    legend(ll2,{'data','model'},'location','NorthWest');
end

end



set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 39 (2+6.5*(length(S)+1))])
print('-depsc2','-r1200',['./figures/' filename]);

