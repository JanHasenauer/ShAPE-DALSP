clear all;
close all;
clc;

%% LOAD DATA
load('./project/results/model_2__alpha_a__beta__a0_delta__grad_on__PL.mat'); parameters_P = parameters;
load('./project/results/model_2__alpha_a__beta__a0_delta__grad_on__MCMC.mat');

%% VISUALIZATION OF PARAMETER DISTRIBUTION
% n_alpha
i = 2;
x = linspace(0,100,25);

figure;
N = histc(10.^parameters.S.par(i,:),x); 
bar(0.5*(x(1:end-1)+x(2:end)),N(1:end-1)/max(N),1,'FaceColor',0.7*[1,1,1]); hold on;
plot(10.^parameters_P.P(i).par(i,:),parameters_P.P(i).R,'-r','linewidth',2);
plot(10.^parameters_P.MS.par(i,1),1,'or','linewidth',2);
plot(x([1,end]),[1,1]*exp(-chi2inv(0.95,1)/2),'--k','linewidth',1);
xlabel('n_\alpha [-]','fontsize',12);
ylabel('posterior','fontsize',12);

text(16,0.40,'statistical','fontsize',10,'HorizontalAlignment','center');
text(16,0.25,'threshold','fontsize',10,'HorizontalAlignment','center');

set(gca,'Fontsize',12,'xlim',x([1,end]),'xtick',0:20:100,'ylim',[0,1.05],'ytick',[]);
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 7.1 4.6]);
print('-depsc2','-r1200',['./figures/fig5_n_alpha']);

% k_alpha
i = 3;
x = linspace(1.45,1.55,30);

figure;
N = histc(10.^parameters.S.par(i,:),x); 
bar(0.5*(x(1:end-1)+x(2:end)),N(1:end-1)/max(N),1,'FaceColor',0.7*[1,1,1]); hold on;
plot(10.^parameters_P.P(i).par(i,:),parameters_P.P(i).R,'-r','linewidth',2);
plot(10.^parameters_P.MS.par(i,1),1,'or','linewidth',2);
plot(x([1,end]),[1,1]*exp(-chi2inv(0.95,1)/2),'--k','linewidth',1);
xlabel('k_\alpha [1/d]','fontsize',12);
ylabel('posterior','fontsize',12);

set(gca,'Fontsize',12,'xlim',x([1,end]),'xtick',[1.46,1.48,1.5,1.52,1.54],'ylim',[0,1.05],'ytick',[]);
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 7.1 4.6]);
print('-depsc2','-r1200',['./figures/fig5_k_alpha']);


% K_alpha
i = 1;
x = linspace(0.4,0.42,30);

figure;
N = histc(10.^parameters.S.par(i,:),x); 
bar(0.5*(x(1:end-1)+x(2:end)),N(1:end-1)/max(N),1,'FaceColor',0.7*[1,1,1]); hold on;
plot(10.^parameters_P.P(i).par(i,:),parameters_P.P(i).R,'-r','linewidth',2);
plot(10.^parameters_P.MS.par(i,1),1,'or','linewidth',2);
plot(x([1,end]),[1,1]*exp(-chi2inv(0.95,1)/2),'--k','linewidth',1);
xlabel('K_\alpha [d]','fontsize',12);
ylabel('posterior','fontsize',12);

set(gca,'Fontsize',12,'xlim',x([1,end]),'xtick',[0.4:0.005:0.42],'ylim',[0,1.05],'ytick',[]);
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 7.1 4.6]);
print('-depsc2','-r1200',['./figures/fig5_K_alpha']);


% k_beta
i = 4;
x = linspace(0.185,0.22,30);

figure;
N = histc(10.^parameters.S.par(i,:),x); 
bar(0.5*(x(1:end-1)+x(2:end)),N(1:end-1)/max(N),1,'FaceColor',0.7*[1,1,1]); hold on;
plot(10.^parameters_P.P(i).par(i,:),parameters_P.P(i).R,'-r','linewidth',2);
plot(10.^parameters_P.MS.par(i,1),1,'or','linewidth',2);
plot(x([1,end]),[1,1]*exp(-chi2inv(0.95,1)/2),'--k','linewidth',1);
xlabel('k_\beta [1/d]','fontsize',12);
ylabel('posterior','fontsize',12);

set(gca,'Fontsize',12,'xlim',x([1,end]),'xtick',[0.19:0.01:0.22],'ylim',[0,1.05],'ytick',[]);
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 7.1 4.6]);
print('-depsc2','-r1200',['./figures/fig5_k_beta']);

%% SIMULATION OF DYNAMICAL SYSTEM DISTRIBUTION

J = 1000;

% Initialization
t_sim = linspace(D.t(1),D.t(end),100)';
I = randperm(size(parameters.S.par,2),J);

pa = nan(length(a_sim),J);
Ni = nan(length(t_sim),8,J);
N = nan(length(t_sim),J);

% Loop: Samples
for i = 1:J
    disp([num2str(100*i/J) '%']);
    % Parameter
    theta_i = parameters.S.par(:,I(i));
    % Simulation
    S = CPsimulateDALSP(M,theta_i,t_sim,a_sim,...
                        [D(1).bins(1,1); D(1).bins(:,2)]',...
                         options_logL.simulation);
    % Save 
    Na_i = sum(S.n_ai(end,:,:),3)';
    I_Na_i = sum(0.5*(Na_i(1:end-1)+Na_i(2:end))*(a_sim(2)-a_sim(1)));
    pa(:,i) = Na_i/I_Na_i;
    Ni(:,:,i) = S.N_i;
    N(:,i) = S.N;
end

% Percentile calculation
pa_perc = prctile(pa',[0.5,50,99.5])';
Ni_perc = nan(length(t_sim),8,3);
N_perc = nan(length(t_sim),3);
for k = 1:length(t_sim)
    N_perc(k,:) = prctile(N(k,:),[0.5,50,99.5]);
    for i = 1:8
        Ni_perc(k,i,:) = prctile(Ni(k,i,:),[0.5,50,99.5]);
    end
end

%
t_sim_plot = t_sim + 1;

%% VISUALIZATION OF SUBPOPULATION SIZE
figure;

Colmap = 0.9*min(colormap(jet(8))+0.2,1);
Colmap_2 = 1.0*min(colormap(jet(8))+0.2,1);

% full figure
lw = 2;

fill([t_sim_plot(1:end);t_sim_plot(end:-1:1)],...
     [N_perc(1:end,1);N_perc(end:-1:1,3)],'k',...
     'FaceColor',0.7*[1,1,1],'EdgeColor',0.7*[1,1,1]); hold on;
legH_Ni(9) = plot(t_sim_plot,N_perc(:,2),'k-','linewidth',lw);
for i = 1:8
    fill([t_sim_plot(1:end);t_sim_plot(end:-1:1)],...
         [Ni_perc(1:end,i,1);Ni_perc(end:-1:1,i,3)],'k',...
         'FaceColor',Colmap_2(i,:),'EdgeColor',Colmap_2(i,:)); hold on;
    legH_Ni(i) = plot(t_sim_plot,Ni_perc(:,i,2),'-','linewidth',lw,'color',Colmap(i,:));
end
xlabel('time [d]','fontsize',12);
ylabel('cell count [-]','fontsize',12);

legend(legH_Ni,{'no division','1 division','2 divisions',...
                '3 divisions','4 divisions','5 divisions',...
                '6 divisions','7 divisions','all'},...
                'Location','NorthEastOutside');
legend('boxoff');

%
xl = 4.8;
xh = 4.9;
yl = 1.1e4;
yh = 1.6e4;
plot([xl,xh,xh,xl,xl],[yl,yl,yh,yh,yl],'k-','linewidth',2);

plot([0.99*xl,1.7],[1.02*yl,3e4],'--k','linewidth',1);
plot([0.99*xh,3.5],[1.02*yh,5.4e4],'--k','linewidth',1);

set(gca,'Fontsize',12,'xlim',D.t([1,end])+1,'xtick',D.t+1,'ylim',[0,6e4]);

% zoom
lw = 4;

axes('position',[0.21,0.56,0.25,0.3]);


fill([t_sim_plot(1:end);t_sim_plot(end:-1:1)],...
     [N_perc(1:end,1);N_perc(end:-1:1,3)],'k',...
     'FaceColor',0.7*[1,1,1],'EdgeColor',0.7*[1,1,1]); hold on;
legH_Ni(9) = plot(t_sim_plot,N_perc(:,2),'k-','linewidth',lw);
for i = 1:8
    fill([t_sim_plot(1:end);t_sim_plot(end:-1:1)],...
         [Ni_perc(1:end,i,1);Ni_perc(end:-1:1,i,3)],'k',...
         'FaceColor',Colmap_2(i,:),'EdgeColor',Colmap_2(i,:)); hold on;
    legH_Ni(i) = plot(t_sim_plot,Ni_perc(:,i,2),'-','linewidth',lw,'color',Colmap(i,:));
end
xlabel('time [d]','fontsize',9);
ylabel('cell count [-]','fontsize',9);
text(0.5*(xl+xh),yl+0.9*(yh-yl),'zoom','fontsize',12,'HorizontalAlignment','center');

set(gca,'Fontsize',9,'xlim',[xl,xh],'xtick',xl:0.02:xh,'ylim',[yl,yh]);

% print
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 18 10]);
print('-depsc2','-r1200',['./figures/fig5_Ni']);

%% VISUALIZATION OF AGE STRUCTURE

figure;

Colmap = [0,0,1];
Colmap_2 = [0.6,0.6,1];
Colmap_3 = [102,255,204]/255;

K_alpha = 10.^parameters.MS.par(1,1);

% full figure
lw = 2;

plot(K_alpha*[1,1],1.7*[0,1],'-','linewidth',lw,'color',Colmap_3); hold on;

fill([a_sim(1:end)';a_sim(end:-1:1)'],...
     [pa_perc(1:end,1);pa_perc(end:-1:1,3)],'k',...
     'FaceColor',Colmap_2,'EdgeColor',Colmap_2); hold on;
plot(a_sim,pa_perc(:,2),'-','linewidth',lw,'color',Colmap);

xlabel('cell age [d]','fontsize',12);
ylabel('age distribution, p(a) [-]','fontsize',12);

set(gca,'xtick',0:0.5:2,'xlim',[0,2],'ylim',[0,1.7],'Fontsize',12);

%
xl = 0.35;
xh = 0.45;
yl = 0.99;
yh = 1.13;
plot([xl,xh,xh,xl,xl],[yl,yl,yh,yh,yl],'k-','linewidth',2);

plot([1.02*xl,1],[1.02*yh,1.6],'--k','linewidth',1);
plot([1.02*xl,1],[0.98*yl,0.72],'--k','linewidth',1);

% zoom
lw = 4;

axes('position',[0.55,0.45,0.33,0.43]);

plot(K_alpha*[1,1],1.7*[0,1],'-','linewidth',lw,'color',Colmap_3); hold on;

fill([a_sim(1:end)';a_sim(end:-1:1)'],...
     [pa_perc(1:end,1);pa_perc(end:-1:1,3)],'k',...
     'FaceColor',Colmap_2,'EdgeColor',Colmap_2); hold on;
plot(a_sim,pa_perc(:,2),'-','linewidth',lw,'color',Colmap);

xlabel('cell age [d]','fontsize',9);
ylabel('age distribution, p(a) [-]','fontsize',9);

set(gca,'Fontsize',9,'xlim',[xl,xh],'xtick',xl:0.02:xh,'ylim',[yl,yh]);

% print
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 13 7]);
print('-depsc2','-r1200',['./figures/fig5_pa']);

