clear all;
close all;
clc;

%%
%load('project/results/model_2__alpha_a__beta__a0_delta__grad_on');
%load('project/results/model_12__alpha_ia__beta_i__a0_delta__grad_on');
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on');

a_max = 4;
pdf.max = [6,6,6,6,6,6];%[];%6;
rate_max = [];%20;
S = 6;%M.S;

%%
a_sim = options_logL.simulation.a_sim;
t_sim = 0;
da = a_sim(2)-a_sim(1);
dim_a_sim = length(a_sim);

%% INITIALIZATION
theta = parameters.MS.MAP.par;

% Bounds
if (~isempty(rate_max)) && (length(rate_max) < S)
    rate_max = rate_max(1) * ones(S,1);
end
if (~isempty(pdf.max)) && (length(pdf.max) < S)
    pdf.max = pdf.max(1) * ones(S,1);
end

% Densities
pdf.A = zeros(length(a_sim),S);
pdf.B = zeros(length(a_sim),S);
pdf.C = zeros(length(a_sim),S);
perc.A = zeros(length(a_sim),S);

% Colors
lblue  = [0.2 0.2 1.0];
lred   = [1.0 0.2 0.2];
lmix   = [0.0 0.6 0.0];
lred_area = [1.000 0.9 0.9];
lblue_area= [0.9 0.9 1.000];
lmix_area = 1 - (1-lblue_area + 1-lred_area);

% Figures
fh1 = figure;
fh2 = figure;

% Loop: subpopulations
for i = 1:S
    figure(fh1);
    
    alpha_i = M.alpha_fun{i}(0,a_sim,theta);
    beta_i  = M.beta_fun{i}(0,a_sim,theta);
    
    %% PLOT DIVISION RATE
    subplot(S,3,3*i-2);
    hl1(1) = plot(a_sim,alpha_i,'LineStyle','-','color',lred,'LineWidth',1.5); hold on;
	hl1(2) = plot(a_sim,beta_i,'LineStyle','-','color',lblue,'LineWidth',1.5); hold on;
	hl1(3) = plot(a_sim,alpha_i+beta_i,'LineStyle','-','color',lmix,'LineWidth',1.5); hold on;
    % Label
    ylabel('rate');
    xlim([0,a_max]);
    if ~isempty(rate_max)
        ylim([0,rate_max(i)]);
    end
    if i == 1
        legend(hl1,{'\alpha_i(a)','\beta_i(a)','\alpha_i(a) + \beta_i(a)'});
    end
    if i ~= S
        set(gca,'xtick',[]);
    else
        xlabel('cell age');
    end
        
    %% COMPUTATION OF PROBABILITY DENSITIES OF INDIVIDUAL EVENTS
    I_a  = [0,0.5*da*cumsum(alpha_i(1:end-1) + alpha_i(2:end))];
    I_b  = [0,0.5*da*cumsum(beta_i(1:end-1)  + beta_i(2:end)) ];
    I_ab = I_a + I_b;
    expI_a   = exp(-I_a);
    expI_b   = exp(-I_b);
    expI_ab  = exp(-I_ab);
    pdf.A(:,i) = alpha_i'.*expI_a';
    pdf.B(:,i) = beta_i'.*expI_b';
    
%     % Loop: cell age
%     for j = 1:length(a)
%         % Division:
%         pdf.A(j,i) = M.alpha_fun{i}(0,a(j),theta)*...
%                           exp(-quad(@(a) M.alpha_fun{i}(0,a,theta),0,a(j)));
%         % Death:
%         pdf.B(j,i) = M.alpha_fun{i}(0,a(j),theta)*...
%                           exp(-quad(@(a) M.beta_fun{i}(0,a,theta),0,a(j)));
%     end
    
    %% PLOT OF PROBABILITY DENSITY OF INDIVIDUAL EVENTS
    % Plot
    subplot(S,3,3*i-1);
    area(a_sim,pdf.A(:,i),'FaceColor',lred_area); hold on;
    area(a_sim,pdf.B(:,i),'FaceColor',lblue_area); hold on;
    area(a_sim,min(pdf.A(:,i),pdf.B(:,i)),'FaceColor',lmix_area); hold on;
    hl2(1) = plot(a_sim,pdf.A(:,i),'LineStyle','-','color',lred,'LineWidth',1.5); hold on;
	hl2(2) = plot(a_sim,pdf.B(:,i),'LineStyle','-','color',lblue,'LineWidth',1.5); hold on;
    % Label
    ylabel('density');
    xlim([0,a_max]);
    if ~isempty(pdf.max)
        ylim([0,pdf.max(i)]);
    end
    if i == 1
        legend(hl2,{'p(a,A|\alpha_i,\beta_i=0)','p(a,B|\alpha_i=0,\beta_i)'});
    end
    if i ~= S
        set(gca,'xtick',[]);
    else
        xlabel('cell age');
    end
    
    %% COMPUTATION OF PROBABILITY DENSITIES OF COMBINED EVENTS
    % Loop: cell age
    for j = 1:length(a_sim)
        % Any event:
        if j == 1
            pdf.C(j,i) =   pdf.A(j,i) ...
                         + pdf.B(j,i) ;
        else
            pdf.C(j,i) =   pdf.A(j,i)*( 1 - 0.5*da*sum(pdf.B(1:j-1,i)+pdf.B(2:j,i)) ) ...
                         + pdf.B(j,i)*( 1 - 0.5*da*sum(pdf.A(1:j-1,i)+pdf.A(2:j,i)) ) ;
            pdf.AB(j,i) = pdf.A(j,i)*( 1 - 0.5*da*sum(pdf.B(1:j-1,i)+pdf.B(2:j,i)) );
            pdf.BA(j,i) = pdf.B(j,i)*( 1 - 0.5*da*sum(pdf.A(1:j-1,i)+pdf.A(2:j,i)) );
        end
    end
    % Division:
%    pdf.AB = pdf.A./(pdf.A + pdf.B).*pdf.C;
%    pdf.BA = pdf.B./(pdf.A + pdf.B).*pdf.C;
%    pdf.AB(:,i) = alpha_i'.*expI_ab';
%    pdf.BA(:,i) = beta_i' .*expI_ab';
    
    %% PLOT OF PROBABILITY DENSITY OF COMBINED EVENTS
    % Plot
    subplot(S,3,3*i-0);
    area(a_sim,pdf.AB(:,i),'FaceColor',lred_area); hold on;
    area(a_sim,pdf.BA(:,i),'FaceColor',lblue_area); hold on;
    area(a_sim,min(pdf.AB(:,i),pdf.BA(:,i)),'FaceColor',lmix_area); hold on;
    hl3(1) = plot(a_sim,pdf.AB(:,i),'LineStyle','-','color',lred,'LineWidth',1.5); hold on;
	hl3(2) = plot(a_sim,pdf.BA(:,i),'LineStyle','-','color',lblue,'LineWidth',1.5); hold on;
    hl3(3) = plot(a_sim,pdf.C(:,i),'LineStyle','-','color',lmix,'LineWidth',1.5); hold on;
    % Label
    ylabel('density');
    xlim([0,a_max]);
    if ~isempty(pdf.max)
        ylim([0,pdf.max(i)]);
    end
    if i == 1
        legend(hl3,{'p(a,A|\alpha_i,\beta_i)','p(a,B|\alpha_i,\beta_i)','p(a,C|\alpha_i,\beta_i)'});
    end
    if i ~= S
        set(gca,'xtick',[]);
    else
        xlabel('cell age');
    end  
    
    %% COMPUTATION OF PERCENTAGE OF DIVIDING AND DYING CELLS
    P.A(i) = 0.5*da*sum(pdf.AB(1:end-1,i) + pdf.AB(2:end,i));    
    for j = 2:length(a_sim)
        perc.A(j,i) = sum(0.5*da*(pdf.AB(1:j-1,i) + pdf.AB(2:j,i)))/P.A(i);
    end
    
    P.B(i) = 0.5*da*sum(pdf.BA(1:end-1,i) + pdf.BA(2:end,i));
    for j = 2:length(a_sim)
        perc.B(j,i) = sum(0.5*da*(pdf.BA(1:j-1,i) + pdf.BA(2:j,i)))/P.B(i);
    end

    %% COMPUTATION MEAN INTER-DIVISION TIME
    MeanA(i) = sum(0.5*da*(a_sim(1:end-1)'.*pdf.AB(1:end-1,i) + a_sim(2:end)'.*pdf.AB(2:end,i)))/P.A(i);
    SigmaA(i) = sqrt(sum(0.5*da*(((a_sim(1:end-1)'-MeanA(i)).^2).*pdf.AB(1:end-1,i) + ((a_sim(2:end)'-MeanA(i)).^2).*pdf.AB(2:end,i))))/P.A(i);

    P.A2(i) = sum(0.5*da*(pdf.A(1:end-1,i) + pdf.A(2:end,i)));    
    MeanA2(i) = sum(0.5*da*(a_sim(1:end-1)'.*pdf.A(1:end-1,i) + a_sim(2:end)'.*pdf.A(2:end,i)))/P.A2(i);
    SigmaA3(i) = sqrt(sum(0.5*da*(((a_sim(1:end-1)'-MeanA(i)).^2).*pdf.A(1:end-1,i) + ((a_sim(2:end)'-MeanA(i)).^2).*pdf.A(2:end,i))))/P.A2(i);

end

%% PLOT SUBPOPULATION STATISTICS
figure(fh2);

% Inter-division time
subplot(2,1,1);
lh4(1) = plot(0:S-1,MeanA,'bo','markersize',8,'linewidth',1.5); hold on;
plot(0:S-1,MeanA,'b--','markersize',8,'linewidth',1); hold on;
for i = 1:S
    ind_05 = min(find(perc.A(:,i) > 0.01))-1;
    ind_50 = min(find(perc.A(:,i) > 0.50))-1;
    ind_95 = max(find(perc.A(:,i) < 0.99))+1;
    % plot(i-1,a(ind_50),'bx','markersize',8,'linewidth',1.5); hold on;
    lh4(2) = plot((i-1)*[1,1],a_sim([ind_05,ind_95]),'b-','linewidth',1.5);
    % plot((i-1)*[1,1],MeanA(i)+SigmaA(i)*[-1,+1],'b-','linewidth',1.5);
end
ylabel('inter-division time');
legend(lh4,{'expected value','5th - 95th percentile'})
xlim([-0.2,S-0.8]);
ylim([0,a_max]);
set(gca,'xtick',[0:S]);

% Size of dividing/dying subpopulation
subplot(2,1,2);
lh5(1) = plot(0:S-1,100*P.A,'bo','markersize',8,'linewidth',1.5); hold on;  
plot(0:S-1,100*P.A,'b--','markersize',8,'linewidth',1); hold on;
lh5(2) = plot(0:S-1,100*P.B,'ro','markersize',8,'linewidth',1.5); hold on;
plot(0:S-1,100*P.B,'r--','markersize',8,'linewidth',1); hold on;
xlabel('subpopulation i');
ylabel('cells [%]');
legend(lh5,{'dividing cells','dying cells'},'Location','Best');
xlim([-0.2,S-0.8]);
ylim([0,100]);
set(gca,'xtick',[0:S]);

