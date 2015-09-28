clear all;
close all;
clc;

%%
markersize = 3;
sig = -(-chi2inv(0.95,1)/2-1);
fontsize = 12;

Imax = 250;

Col = [255,128,0;...
       255,204,102;...
       0,128,255;... % 102,204,255;...
       100,100,100]/255;

Coltrans = [255,200,200;...
            255,240,190;...
            200,200,255;... % 102,204,255;...
            100,100,100]/255;

%%
name = {'random sampling',...
        ['stochastic optimisation' char(10) '(particle-swarm)'],...
        ['deterministic optimisation' char(10) '(finite differences)'],...
        ['deterministic optimisation' char(10) '(sensitivity equations)']};

%% Model 1
% Load data
load('project/results/model_1__alpha__beta__a0_delta__grad_on');                      params{1}  = parameters; II(1)  = 1;
load('project/results/model_1__alpha__beta__a0_delta__grad_off');                     params{2}  = parameters; II(2)  = 2;
load('project/results/model_1__alpha__beta__a0_delta__PS__Size_120__MaxObj_12000');   params{3}  = parameters; II(3)  = 3;

% random sampling
load('project/results/model_1__alpha__beta__a0_delta__grad_off');
parameters.MS.logPost = sort(parameters.MS.logPost0,1,'descend');
parameters.MS.logPost = [parameters.MS.logPost(~isnan(parameters.MS.logPost));parameters.MS.logPost(isnan(parameters.MS.logPost))];
params{4}  = parameters; II(4)  = 4;

% Objective function
figure;
axes('Position',[0.1,0.47,0.9,0.48]);
%axes('Position',[0.1,0.47,0.867,0.48]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

plot([1-10,Imax+10],(sig+1)*[1,1],'k--');
text(165,9,'statistical threshold','fontsize',fontsize-2);

set(gca,'xtick',[0,50,100,150,200,250],'xlim',[1-10,Imax+10]);
set(gca,'ytick',[1,1e2,1e4,1e6,1e8,1e10,1e12],'yscale','log','ylim',[0.5,1e12]);

xlabel('sorted optimiser runs','fontsize',fontsize)
ylabel('negative log-likelihood + constant','fontsize',fontsize);

legend(lh(end:-1:1),name,'location','NorthEastOutside','fontsize',fontsize);
legend('boxoff')

% 
axes('Position',[0.15,0.76,0.28,0.17]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

set(gca,'xtick',[0,50,100,150],'xlim',[-4,154],'fontsize',fontsize-4);
set(gca,'yscale','log','ylim',[0.99,1.1],'fontsize',fontsize-4);

xlabel('sorted optimiser runs','fontsize',fontsize-4)
text(75,1.085,'Zoom','fontsize',fontsize-4)

% Computation time
axes('Position',[0.3383,0.12,0.28,0.22]);

options_boxplot.xscale = 'log';
options_boxplot.line_color = Col(end-1:-1:1,:);
options_boxplot.fill_color = Coltrans(end-1:-1:1,:);
options_boxplot.marker_size = markersize;
boxplot_horizontal({params{3}.MS.t_cpu/60^2,params{2}.MS.t_cpu/60^2,params{1}.MS.t_cpu/60^2},{'','',''},options_boxplot);

xl = [3*1e-2,20e1];
x_text = 10.^(log10(xl(1))-(log10(xl(2))-log10(xl(1)))*0.0587);
text(x_text,0,name{4},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,1,name{3},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,2,name{2},'HorizontalAlignment','right','fontsize',fontsize);

set(gca,'xscale','log','xtick',[0.01,0.1,1,10,100],'xlim',xl,'ylim',[-0.5,2.5]);
xlabel({'computation time','for individual start [hr]'});

% Computation time
axes('Position',[0.3383+0.32,0.12,0.28,0.22]);

logPost_opt = max([max(params{1}.MS.logPost),...
                   max(params{2}.MS.logPost),...
                   max(params{3}.MS.logPost)]);
               
for i = 1:3
    n(i) = sum((logPost_opt - params{i}.MS.logPost) < sig);
    e(i) = nansum(params{i}.MS.t_cpu/60^2)/n(i);
end
    
for i = 1:3
    plot(e(i),i,'o','color',Col(i,:),'MarkerFaceColor',Coltrans(i,:)); hold on;
end

set(gca,'xscale','log','xtick',[1,10,100],'xlim',[0.4e0,0.3e2]);
set(gca,'ytick',[],'ylim',[0.5,3.5]);
xlabel({'average computation time','per converged start [hr]'});

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 17 17])
print('-depsc2','-r1200',['./figures/figS2a']);

%% Model 2
% Load data
load('project/results/model_2__alpha_a__beta__a0_delta__grad_on');                      params{1}  = parameters; II(1)  = 1;
load('project/results/model_2__alpha_a__beta__a0_delta__grad_off');                     params{2}  = parameters; II(2)  = 2;
load('project/results/model_2__alpha_a__beta__a0_delta__PS__Size_140__MaxObj_14000');   params{3}  = parameters; II(3)  = 3;

% random sampling
load('project/results/model_2__alpha_a__beta__a0_delta__grad_off');
parameters.MS.logPost = sort(parameters.MS.logPost0,1,'descend');
parameters.MS.logPost = [parameters.MS.logPost(~isnan(parameters.MS.logPost));parameters.MS.logPost(isnan(parameters.MS.logPost))];
params{4}  = parameters; II(4)  = 4;

% Objective function
figure;
axes('Position',[0.1,0.47,0.9,0.48]);
%axes('Position',[0.1,0.47,0.867,0.48]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

plot([1-10,Imax+10],(sig+1)*[1,1],'k--');
text(165,9,'statistical threshold','fontsize',fontsize-2);

set(gca,'xtick',[0,50,100,150,200,250],'xlim',[1-10,Imax+10]);
set(gca,'ytick',[1,1e2,1e4,1e6,1e8,1e10,1e12],'yscale','log','ylim',[0.5,1e12]);

xlabel('sorted optimiser runs','fontsize',fontsize)
ylabel('negative log-likelihood + constant','fontsize',fontsize);

legend(lh(end:-1:1),name,'location','NorthEastOutside','fontsize',fontsize);
legend('boxoff')

% 
axes('Position',[0.15,0.76,0.18,0.17]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

set(gca,'xtick',[0,20,40,60],'xlim',[-4,72],'fontsize',fontsize-4);
set(gca,'yscale','log','ylim',[0.99,1.1],'fontsize',fontsize-4);

xlabel('sorted optimiser runs','fontsize',fontsize-4)
text(30,1.085,'Zoom','fontsize',fontsize-4)

% Computation time
axes('Position',[0.3383,0.12,0.28,0.22]);

options_boxplot.xscale = 'log';
options_boxplot.line_color = Col(end-1:-1:1,:);
options_boxplot.fill_color = Coltrans(end-1:-1:1,:);
options_boxplot.marker_size = markersize;
boxplot_horizontal({params{3}.MS.t_cpu/60^2,params{2}.MS.t_cpu/60^2,params{1}.MS.t_cpu/60^2},{'','',''},options_boxplot);

xl = [0.5*1e-2,3e1];
x_text = 10.^(log10(xl(1))-(log10(xl(2))-log10(xl(1)))*0.0587);
text(x_text,0,name{4},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,1,name{3},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,2,name{2},'HorizontalAlignment','right','fontsize',fontsize);

set(gca,'xscale','log','xtick',[0.01,0.1,1,10,100],'xlim',xl,'ylim',[-0.5,2.5]);
xlabel({'computation time','for individual start [hr]'});

% Computation time
axes('Position',[0.3383+0.32,0.12,0.28,0.22]);

logPost_opt = max([max(params{1}.MS.logPost),...
                   max(params{2}.MS.logPost),...
                   max(params{3}.MS.logPost)]);
               
for i = 1:3
    n(i) = sum((logPost_opt - params{i}.MS.logPost) < sig);
    e(i) = nansum(params{i}.MS.t_cpu/60^2)/n(i);
end
    
for i = 1:3
    plot(e(i),i,'o','color',Col(i,:),'MarkerFaceColor',Coltrans(i,:)); hold on;
end

set(gca,'xscale','log','xtick',[1,10,100],'xlim',[0.7*1e0,2e2]);
set(gca,'ytick',[],'ylim',[0.5,3.5]);
xlabel({'average computation time','per converged start [hr]'});

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 17 17])
print('-depsc2','-r1200',['./figures/fig2']);

%% Model 16
% Load data
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on');                      params{1}  = parameters; II(1)  = 1;
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_off');                     params{2}  = parameters; II(2)  = 2;
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__PS__Size_290__MaxObj_29000');   params{3}  = parameters; II(3)  = 3;

% random sampling
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_off');
parameters.MS.logPost = sort(parameters.MS.logPost0,1,'descend');
parameters.MS.logPost = [parameters.MS.logPost(~isnan(parameters.MS.logPost));parameters.MS.logPost(isnan(parameters.MS.logPost))];
params{4}  = parameters; II(4)  = 4;

% Objective function
figure;
axes('Position',[0.1,0.47,0.9,0.48]);
%axes('Position',[0.1,0.47,0.867,0.48]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

plot([1-10,Imax+10],(sig+1)*[1,1],'k--');
text(165,9,'statistical threshold','fontsize',fontsize-2);

set(gca,'xtick',[0,50,100,150,200,250],'xlim',[1-10,Imax+10]);
set(gca,'ytick',[1,1e2,1e4,1e6,1e8,1e10,1e12],'yscale','log','ylim',[0.5,1e12]);

xlabel('sorted optimiser runs','fontsize',fontsize)
ylabel('negative log-likelihood + constant','fontsize',fontsize);

legend(lh(end:-1:1),name,'location','NorthEastOutside','fontsize',fontsize);
legend('boxoff')

% 
axes('Position',[0.15,0.76,0.13,0.17]);

for i = 1:4
    lh(i) = plot(1:Imax,1+params{1}.MS.logPost(1)-params{i}.MS.logPost,'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
end

for j = 1:Imax
    for i = mod([1:4]+j+2,4)+1
        plot(j,1+params{1}.MS.logPost(1)-params{i}.MS.logPost(j),'-o','markersize',markersize,'color',Col(i,:),'MarkerFaceColor',[1,1,1]); hold on;
    end
end

set(gca,'xtick',[0,5,10,15],'xlim',[-2,19],'fontsize',fontsize-4);
set(gca,'yscale','log','ylim',[0.99,1.1],'fontsize',fontsize-4);

xlabel('sorted optimiser runs','fontsize',fontsize-4)
text(7,1.085,'Zoom','fontsize',fontsize-4)

% Computation time
axes('Position',[0.3383,0.12,0.28,0.22]);

options_boxplot.xscale = 'log';
options_boxplot.line_color = Col(end-1:-1:1,:);
options_boxplot.fill_color = Coltrans(end-1:-1:1,:);
options_boxplot.marker_size = markersize;
boxplot_horizontal({params{3}.MS.t_cpu/60^2,params{2}.MS.t_cpu/60^2,params{1}.MS.t_cpu/60^2},{'','',''},options_boxplot);

xl = [0.5*1e-2,0.8e2];
x_text = 10.^(log10(xl(1))-(log10(xl(2))-log10(xl(1)))*0.0587);
text(x_text,0,name{4},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,1,name{3},'HorizontalAlignment','right','fontsize',fontsize);
text(x_text,2,name{2},'HorizontalAlignment','right','fontsize',fontsize);

set(gca,'xscale','log','xtick',[0.01,0.1,1,10,100],'xlim',xl,'ylim',[-0.5,2.5]);
xlabel({'computation time','for individual start [hr]'});

% Computation time
axes('Position',[0.3383+0.32,0.12,0.28,0.22]);

logPost_opt = max([max(params{1}.MS.logPost),...
                   max(params{2}.MS.logPost),...
                   max(params{3}.MS.logPost)]);
               
for i = 1:3
    n(i) = sum((logPost_opt - params{i}.MS.logPost) < sig);
    e(i) = nansum(params{i}.MS.t_cpu/60^2)/n(i);
end
    
for i = 1:3
    plot(e(i),i,'o','color',Col(i,:),'MarkerFaceColor',Coltrans(i,:)); hold on;
end

set(gca,'xscale','log','xtick',[1,10,100,1000],'xlim',[1e1,1e3]);
set(gca,'ytick',[],'ylim',[0.5,3.5]);
xlabel({'average computation time','per converged start [hr]'});

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 17 17])
print('-depsc2','-r1200',['./figures/figS2b']);

