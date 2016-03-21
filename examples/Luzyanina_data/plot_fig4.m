clear all;
close all;
clc;

%%
name = {'model_1__alpha__beta__a0_delta',...
        'model_2__alpha_a__beta__a0_delta',...
        'model_3__alpha_i__beta__a0_delta',...
        'model_4__alpha_ia__beta__a0_delta',...
        'model_5__alpha__beta_a__a0_delta',...
        'model_6__alpha_a__beta_a__a0_delta',...
        'model_7__alpha_i__beta_a__a0_delta',...
        'model_8__alpha_ia__beta_a__a0_delta',...
        'model_9__alpha__beta_i__a0_delta',...
        'model_10__alpha_a__beta_i__a0_delta',...
        'model_11__alpha_i__beta_i__a0_delta',...
        'model_12__alpha_ia__beta_i__a0_delta',...
        'model_13__alpha__beta_ia__a0_delta',...
        'model_14__alpha_a__beta_ia__a0_delta',...
        'model_15__alpha_i__beta_ia__a0_delta',...
        'model_16__alpha_ia__beta_ia__a0_delta'};

Dep = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
       0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1
       0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1
       0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];

% Assign data
n_theta_vec = nan(length(name),1);
logL_vec = nan(length(name),1);
AIC_vec = nan(length(name),1);
BIC_vec = nan(length(name),1);
for i = 1:length(name)
    pathname = ['./project/results/' name{i} '__grad_on'];
    load(pathname);
    n_theta_vec(i) = parameters.number;
    logL_vec(i) = parameters.MS.logPost(1);
    AIC_vec(i) = -2*parameters.MS.logPost(1) + 2*parameters.number;
    BIC_vec(i) = -2*parameters.MS.logPost(1) + parameters.number*log(sum(sum(D.cellcount)));
end


%%
figure;

% Determine index set

x_lim = -[2.45,1.15]*1e4/1e4;
dx = diff(x_lim)/30;
dy = 0*(max(n_theta_vec) - min(n_theta_vec))/20;

CC = [102,204,255;...
      0.8*102,255,0.8*204;...
      1.2*204,255,1.2*102;...
      255,204,102;...
      ]/255;

for i = 1:16
    % Mark
    if Dep(1,i) == 1
        plot(logL_vec(i)/1e4,n_theta_vec(i),'o','markersize',13,'color',CC(1,:),'MarkerFaceColor',CC(1,:)); hold on;
        if Dep(2,i) == 1
            plot(logL_vec(i)/1e4,n_theta_vec(i),'o','markersize',10,'color',CC(2,:),'MarkerFaceColor',CC(2,:)); hold on;
            if Dep(4,i) == 1
                plot(logL_vec(i)/1e4,n_theta_vec(i),'o','markersize',7,'color',CC(3,:),'MarkerFaceColor',CC(3,:)); hold on;
                if Dep(3,i) == 1
                    plot(logL_vec(i)/1e4,n_theta_vec(i),'o','markersize',4,'color',CC(4,:),'MarkerFaceColor',CC(4,:)); hold on;
                end
            end
        end
    else
        plot(logL_vec(i)/1e4,n_theta_vec(i),'o','markersize',10,'linewidth',4,'color',0.7*[1,1,1]); hold on;
    end
    % Name
    text((logL_vec(i)/1e4+dx),n_theta_vec(i)+dy,['{\bf $\mathcal{M}_{' num2str(i,'%d') '}$}'],'Interpreter','latex','fontsize',12)
    %
    axis square;
    xlim(x_lim);
    ylim([10,32]);
end

text(-1.58,11.3,...
    ['$\mathrm{P}_{\mathrm{corr}}$ = ' num2str(corr(n_theta_vec,logL_vec))],'Interpreter','latex','fontsize',12);
%
legend({'$\alpha$ is independent of $a$ ($\alpha \neq$ fnc($a$))\\[2ex]',...
        '$\alpha$ is depend on $a$ ($\alpha =$ fnc($a$))'},...
        'position',[0.2295 0.8 0.3172 0.1],'Interpreter','latex');
%legend('\alpha \neq fnc({\ita})','\alpha = fnc({\ita})','location','NorthWest');
%xlabel('BIC [$10^4$]','Interpreter','latex','fontsize',12);
%ylabel('number of model parameters, n_\theta','Interpreter','latex','fontsize',12);

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 14 13])
print('-depsc2','-r2000',['./figures/fig4a']);
print('-dpdf','-r2000',['./figures/fig4a']);


%%
figure;

% Determine index set
[BIC_sorted,I] = sort(BIC_vec);

x_lim = [0.5,16.5];
y_lim = [2.5,4.6]*1e4/1e4;

subplot('Position',[0.22,0.48,0.76,0.45]);
plot(1:16,BIC_sorted/1e4,'k-o','linewidth',1,'markersize',6,'MarkerFaceColor',[0,0,0]); hold on;
set(gca,'xtick',[]);
ylabel('BIC [$10^4$]','Interpreter','latex','fontsize',12);

xlim(x_lim);
ylim(y_lim);

%
subplot('Position',[0.22,0.13,0.76,0.32]);
J = [1,2,4,3];
for i = 1:16
    for j = 1:4
        if Dep(J(j),I(i)) == 1
            plot(i,4-j,'kx','linewidth',1,'markersize',12); hold on;
        end
    end
    if i == 1
        for j = 1:4
            switch J(j)
                case 1
                    text(-2.8,4-j,'$\alpha$ = fnc($a$)','Interpreter','latex','fontsize',12);
                case 2
                    text(-2.8,4-j,'$\alpha$ = fnc($i$)','Interpreter','latex','fontsize',12);
                case 3
                    text(-2.8,4-j,'$\beta$ = fnc($a$)' ,'Interpreter','latex','fontsize',12);
                case 4
                    text(-2.8,4-j,'$\beta$ = fnc($i$)' ,'Interpreter','latex','fontsize',12);
            end
        end
    end
end
xlim(x_lim);
ylim([-0.5,3.5]);

for j = 1:3
    plot([0.5,16.5],-[0.5,0.5]+j,'k-','linewidth',1);
end
for j = 1:15
    plot(j+[0.5,0.5],[-0.5,3.5],'k-','linewidth',1);
end
set(gca,'xtick',[],'ytick',[]);

for i = 1:16
	text(i,-1.5,['{\bf $\mathcal{M}_{' num2str(I(i),'%d') '}$}'],'Interpreter','latex','fontsize',12,'rotation',90);
end

%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 13.2])
print('-depsc2','-r1200',['./figures/fig4b']);
