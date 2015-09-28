% Plot
for i = 1:5
figure;
area(a,pdf.AB(:,i),'FaceColor',lred_area); hold on;
area(a,pdf.BA(:,i),'FaceColor',lblue_area); hold on;
area(a,min(pdf.AB(:,i),pdf.BA(:,i)),'FaceColor',lmix_area); hold on;
hl3(1) = plot(a,pdf.AB(:,i),'LineStyle','-','color',lred,'LineWidth',1.5); hold on;
hl3(2) = plot(a,pdf.BA(:,i),'LineStyle','-','color',lblue,'LineWidth',1.5); hold on;
hl3(3) = plot(a,pdf.C(:,i),'LineStyle','-','color',lmix,'LineWidth',1.5); hold on;
% Label
xlabel('cell age');
ylabel('density');
xlim([0,a_max]);
if ~isempty(pdf.max)
    ylim([0,pdf.max(i)]);
end
set(gcf, 'Position',         [0 0 7 4] );
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc',['WCSB_DALSP_pdfs_' num2str(i-1)])

end

figure;
set(gcf, 'Position',         [0 0 11 12] );

% Inter-division time
subplot(2,1,1);
lh4(1) = plot(0:S-1,MeanA,'bo','markersize',8,'linewidth',1.5); hold on;
plot(0:S-1,MeanA,'b--','markersize',8,'linewidth',1); hold on;
for i = 1:S
    ind_05 = min(find(perc.A(:,i) > 0.05))-1;
    ind_50 = min(find(perc.A(:,i) > 0.50))-1;
    ind_95 = max(find(perc.A(:,i) < 0.95))+1;
    % plot(i-1,a(ind_50),'bx','markersize',8,'linewidth',1.5); hold on;
    lh4(2) = plot((i-1)*[1,1],a([ind_05,ind_95]),'b-','linewidth',1.5);
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

set(gcf, 'PaperPositionMode', 'auto');
print('-depsc','WCSB_DALSP_statistics')

figure(options_logL.figure_handle)
set(gcf, 'Position',  [0 0 20 10] );
logLikelihood_proliferation(parameters.ml,M,D,options_logL);
xlim([5,3000])
xlabel('measured fluorescence','fontsize',18)
ylabel('cell count','fontsize',18)
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc','WCSB_DALSP_fit')

%%
figure;
set(gcf, 'Position',         [0 0 210 130] );
AA = 0:0.01:4;
k = 0.8;
y1 = k*exp(-k*AA);
y2 = lognpdf(AA,0.1,0.3);

area(AA,y1,'FaceColor',lred_area); hold on;
% area(AA,y2,'FaceColor',lblue_area); hold on;
% area(AA,min(y1,y2),'FaceColor',lmix_area); hold on;
plot(AA,k*exp(-k*AA),'-','color',lred,'linewidth',1.5); hold on;
% plot(AA,y2,'-','color',lblue,'linewidth',1.5); hold on;



xlabel('inter-division time','fontsize',14)
ylabel('density','fontsize',14)
ylim([0,1.5])

set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'xtick',[]); set(gca, 'ytick',[]);

print('-depsc',['WCSB_DALSP_pdfs_illustration_1'])

%%
figure;
set(gcf, 'Position',         [0 0 210 130] );
AA = 0:0.01:4;
k = 0.8;
y1 = k*exp(-k*AA);
y2 = lognpdf(AA,0.1,0.3);

area(AA,y1,'FaceColor',lred_area); hold on;
area(AA,y2,'FaceColor',lblue_area); hold on;
area(AA,min(y1,y2),'FaceColor',lmix_area); hold on;
plot(AA,k*exp(-k*AA),'-','color',lred,'linewidth',1.5); hold on;
plot(AA,y2,'-','color',lblue,'linewidth',1.5); hold on;



xlabel('inter-division time','fontsize',14)
ylabel('density','fontsize',14)
ylim([0,1.5])

set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'xtick',[]); set(gca, 'ytick',[]);

print('-depsc',['WCSB_DALSP_pdfs_illustration_2'])
