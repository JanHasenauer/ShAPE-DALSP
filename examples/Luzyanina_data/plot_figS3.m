clear all;
close all;
clc;

%% LOAD DATA
load('project/results/model_1__alpha__beta__a0_delta__grad_on');        Mo{1}  = M; par{1}  = parameters; II(1)  = 1;
load('project/results/model_2__alpha_a__beta__a0_delta__grad_on');      Mo{2}  = M; par{2}  = parameters; II(2)  = 2;
load('project/results/model_3__alpha_i__beta__a0_delta__grad_on');      Mo{3}  = M; par{3}  = parameters; II(3)  = 3;
load('project/results/model_4__alpha_ia__beta__a0_delta__grad_on');     Mo{4}  = M; par{4}  = parameters; II(4)  = 4;
load('project/results/model_5__alpha__beta_a__a0_delta__grad_on');      Mo{5}  = M; par{5}  = parameters; II(5)  = 5;
load('project/results/model_6__alpha_a__beta_a__a0_delta__grad_on');    Mo{6}  = M; par{6}  = parameters; II(6)  = 6;
load('project/results/model_7__alpha_i__beta_a__a0_delta__grad_on');    Mo{7}  = M; par{7}  = parameters; II(7)  = 7;
load('project/results/model_8__alpha_ia__beta_a__a0_delta__grad_on');   Mo{8}  = M; par{8}  = parameters; II(8)  = 8;
load('project/results/model_9__alpha__beta_i__a0_delta__grad_on');      Mo{9}  = M; par{9}  = parameters; II(9)  = 9;
load('project/results/model_10__alpha_a__beta_i__a0_delta__grad_on');   Mo{10} = M; par{10} = parameters; II(10) = 10;
load('project/results/model_11__alpha_i__beta_i__a0_delta__grad_on');   Mo{11} = M; par{11} = parameters; II(11) = 11;
load('project/results/model_12__alpha_ia__beta_i__a0_delta__grad_on');  Mo{12} = M; par{12} = parameters; II(12) = 12;
load('project/results/model_13__alpha__beta_ia__a0_delta__grad_on');    Mo{13} = M; par{13} = parameters; II(13) = 13;
load('project/results/model_14__alpha_a__beta_ia__a0_delta__grad_on');  Mo{14} = M; par{14} = parameters; II(14) = 14;
load('project/results/model_15__alpha_i__beta_ia__a0_delta__grad_on');  Mo{15} = M; par{15} = parameters; II(15) = 15;
load('project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on'); Mo{16} = M; par{16} = parameters; II(16) = 16;

%% 
figure;

Col = [255,128,0;...
       255,204,102;...
       0,128,255;... % 102,204,255;...
       100,100,100]/255;

getIsp = @(i,j) [(7*(i-1)+2+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)],...
                 (7*(i-1)+3+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)],...
                 (7*(i-1)+4+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)],...
                 (7*(i-1)+5+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)],...
                 (7*(i-1)+6+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)],...
                 (7*(i-1)+7+6)*(4*2+2)+[(4*(j-1)+1)+(j>1):4*j+(j>1)]];

             
Imax = 250;
for i = 1:16
        %
        subplot(8+8*(7+1),4*2+2,getIsp(mod(i-1,8)+1,1+(i>8)));
%        plot(1:Imax,[par{i}.MS.logPost(1:min(Imax,length(par{i}.MS.logPost))),nan(1,Imax-length(par{i}.MS.logPost))],'k-o','markersize',2);
        plot(1:Imax,1+par{16}.MS.logPost(1)-[par{i}.MS.logPost(1:min(Imax,length(par{i}.MS.logPost))),nan(1,Imax-length(par{i}.MS.logPost))],'-o',...
            'color',0*[1,1,1],'markersize',3,'MarkerFaceColor',[1,1,1]); hold on;
        %
        xlim([1,Imax]);
%        ylim([-10e4,-1e4]);
        set(gca,'ytick',[1,1e2,1e4,1e6],'yscale','log','ylim',[0.5,1e6]);
        if i == 8 || i == 16
            set(gca,'xtick',0:50:250);
            xlabel('sorted optimiser runs');
        else
            set(gca,'xtick',0:50:250,'xticklabel',[]);
        end
%        ylabel('log-likelihood');
        ylabel({'neg. log-likelihood',' + constant'});
%        xA = -26*250/100; yA = -0.3e4;
        if i <= 8
            xA = -28*250/100; yA = 3e6;
        else
            xA = -30*250/100; yA = 3e6;
        end
        switch i
            case 1, text(xA,yA,'{\bf A}','fontsize',16);
            case 2, text(xA,yA,'{\bf B}','fontsize',16);
            case 3, text(xA,yA,'{\bf C}','fontsize',16);
            case 4, text(xA,yA,'{\bf D}','fontsize',16);
            case 5, text(xA,yA,'{\bf E}','fontsize',16);
            case 6, text(xA,yA,'{\bf F}','fontsize',16);
            case 7, text(xA,yA,'{\bf G}','fontsize',16);
            case 8, text(xA,yA,'{\bf H}','fontsize',16);
            case 9, text(xA,yA,'{\bf I}','fontsize',16);
            case 10, text(xA,yA,'{\bf J}','fontsize',16);
            case 11, text(xA,yA,'{\bf K}','fontsize',16);
            case 12, text(xA,yA,'{\bf L}','fontsize',16);
            case 13, text(xA,yA,'{\bf M}','fontsize',16);
            case 14, text(xA,yA,'{\bf N}','fontsize',16);
            case 15, text(xA,yA,'{\bf O}','fontsize',16);
            case 16, text(xA,yA,'{\bf P}','fontsize',16);
        end
        text(xA+6*250/100,2.5e6,['{\bf $(\mathcal{M}_{' num2str(i,'%d') '})$}'],'Interpreter','latex','fontsize',16);

        plot([1-3,Imax+3],[1,1],'-','color',Col(1,:)); xlim([1-3,Imax+3]);

%         sig = -(-chi2inv(0.95,1)/2-1);
%         plot([1-3,Imax+3],sig*[1,1],'r--');

        if i == 16
            text(160,4,'optimum for','fontsize',12,'color',Col(1,:));
            text(160+49,3.7,'$\mathcal{M}_{16}$','Interpreter','latex','fontsize',12,'color',Col(1,:));
%             text(165,7,'statistical threshold','fontsize',10,'color',[1,0,0]);
        end
end


% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 41 70])
% print('-depsc2','-r1200',['./figures/figS3']);
