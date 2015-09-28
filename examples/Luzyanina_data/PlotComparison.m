%% Compare estimations

P1 = load('project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P','parameters','options_PL');
P2 = load(['project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P_reopt.mat'],'parameters','options_PL');
P3 = load(['project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P_new.mat'],'parameters','options_PL');
%P3 = load(['EstResults/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P',num2str(k),'.mat'],'parameters','options_PL');
%P3 = load(['project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P',num2str(k),'.mat'],'parameters','options_PL');
% P3 = load('2013-09-26/ProfTemp-08-44.mat','parameters','options');



figure(1);
for k = P1.options_PL.parameter_index
    ind = find(P2.parameters.MS.MAP_list.logPost >= P2.parameters.MS.MAP.logPost-chi2inv(0.95,29)/2);
    subplot(5,6,k)
    hold on;
    plot(P1.parameters.P(k).par(k,:),exp(P1.parameters.P(k).logPost-P1.parameters.MS.MAP.logPost),'b','Linewidth',1.5); 
    plot(P2.parameters.P(k).par(k,:),exp(P2.parameters.P(k).logPost-P2.parameters.MS.MAP.logPost),'+-k','Linewidth',1.5); 
    xl = [min([P2.parameters.P(k).par(k,:),min(P2.parameters.MS.MAP_list.par(k,ind))]),...
          max([P2.parameters.P(k).par(k,:),max(P2.parameters.MS.MAP_list.par(k,ind))])];
    xlim(xl);
    ylim([0 1]);
end

for k = P3.options_PL.parameter_index
    ind = find(P3.parameters.MS.MAP_list.logPost >= P3.parameters.MS.MAP.logPost-chi2inv(0.95,29)/2);
    subplot(5,6,k)
    hold on;
    plot(P3.parameters.P(k).par(k,:),exp(P3.parameters.P(k).logPost-P3.parameters.MS.MAP.logPost),'g','Linewidth',1.5); hold on;
    xl = [min([P3.parameters.P(k).par(k,:),min(P3.parameters.MS.MAP_list.par(k,ind))]),...
          max([P3.parameters.P(k).par(k,:),max(P3.parameters.MS.MAP_list.par(k,ind))])];
    xlim(xl);
    ylim([0 1]);
end
legend('trust-region-reflective TolFun 1e-6','reoptimized', 'trust-region-reflective TolFun 1e-8')

% figure(2)
% for k = ind
%     subplot(5,6,k)
%     plot(P1.parameters.P(k).par(k,:),exp(P1.parameters.P(k).logPost-P1.parameters.MS.MAP.logPost),'b','Linewidth',1.5); hold on;
%     plot(P2.parameters.P(k).par(k,:),exp(P2.parameters.P(k).logPost-P2.parameters.MS.MAP.logPost),'k','Linewidth',1.5); hold on;
%     plot(P3.parameters.P(k).par(k,:),exp(P3.parameters.P(k).logPost-P3.parameters.MS.MAP.logPost),'-og','Linewidth',1.5);
% end
% legend('trust-region-reflective TolFun 1e-6','reoptimized', 'trust-region-reflective TolFun 1e-8')


%% Plot all

P1 = load('project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P','parameters','options_PL');
P2 = load('project/results/model_16__alpha_ia__beta_ia__a0_delta__PL__grad_on_P_reopt.mat','parameters','options_PL');

p = 1e-4;
xmin = (1-p)*P2.parameters.min +    p *P2.parameters.max;
xmax =    p *P2.parameters.min + (1-p)*P2.parameters.max;

figure;
s(1) = 5;
s(2) = 6;

% Loop: Parameter
for i = P2.options_PL.parameter_index
    % Open subplot
    subplot(s(1),s(2),i);
    hold on;

    % Plot profile likelihood
    plot(P1.parameters.P(i).par(i,:),exp(P1.parameters.P(i).logPost - P1.parameters.MS.MAP.logPost),'r-','linewidth',1.5); hold on;
    plot(P2.parameters.P(i).par(i,:),exp(P2.parameters.P(i).logPost - P2.parameters.MS.MAP.logPost),'b-','linewidth',1.5); hold on;
    
    ind = find(P1.parameters.MS.MAP_list.logPost >= P1.parameters.MS.MAP.logPost-chi2inv(0.01,29)/2);
    xl = [min([P1.parameters.P(i).par(i,:),min(P1.parameters.MS.MAP_list.par(i,ind))]),...
          max([P1.parameters.P(i).par(i,:),max(P1.parameters.MS.MAP_list.par(i,ind))])];
      
    xlim(xl); 
    ylim([0,1.1]);
end
