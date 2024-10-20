%                                 _---~~(~~-_.
%                               _{        )   )
% ██████  ███████ ██████      ,   ) -~~- ( ,-' )_
% ██   ██ ██      ██   ██    (  `-,_..`., )-- '_,)
% ██████  █████   ██   ██   ( ` _)  (  -~( -_ `,  }
% ██   ██ ██      ██   ██   (_-  _  ~_-~~~~`,  ,' )
% ██   ██ ███████ ██████      `~ -^(    __;-,((()))
% Richard E. Daws  - 2021           ~~~~ {_ -_(())
%                                          `\  }
%                                            { }   
% Psilodep :- Modularity visualisation
% 
% Requires 
%    third_party folder and subfolders to be in path 


%% Setup

clear; close all

% Boolean to save out figs or not.
saveFigs = false;            


% load in data from each study & plotting variables
load('../data/psilodep1/dat_1.mat');
load('../data/psilodep2/dat_2.mat'); 
load('../data/psilodep2/flex.mat'); 
load('../plotting_vars.mat') 

% Final figure size
fgPos = [1 1 1440 796];
% Data dot size
plt.mrkSze = 50; 
% plot font size
plt.fntSze = 25;            

% Pull the modularity values for each study and session
Q1 = dat_1{:,{'before_rest','after_rest'}};
Q2 = dat_2{:,{'ses_1','ses_2'}};

%% PSILODEP 1 %%

% Define the comparison labels 
cpLbl = {'after_rest:BDI_week_1_post_25mg','after_rest:BDI_month_3_post_25mg','after_rest:BDI_month_6_post_25mg'}; 
% Define the statistic labels
stLbl = {'R','CI','p'};

% Create the output table
T.corr = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
% Add Confidence interval column
T.corr.CI = NaN(height(T.corr), 2);

% Run a correlation for each comparison label 
for ii = 1:numel(cpLbl)
    % Split the string with ":" delim to indentify comparison timepoints
    tmp = split(T.corr.Row(ii),':');
    % Pearson corr - two-tailed
    [R,P] = corr(dat_1{:, tmp(1)}, dat_1{:, tmp(2)});
    % Add results to table with confidence interval
    T.corr{T.corr.Row(ii), :} = [R c3nl_RCI(R,height(dat_1),1) P];
end    

% Add FDR corrected p-values   
[~, ~, ~, Pc]=fdr_bh(T.corr.p, 0.05, 'pdep', 'no');
T.corr.p_fdr = Pc;

% Correlation between the changes in primary endpoint BDI and Q   
T.diff = array2table(NaN(1,numel(stLbl)), 'VariableNames', stLbl,'RowNames',{'Q_BDI_change'});
T.diff.CI = NaN(height(T.diff), 2);

[R,P] = corr(dat_1.after_rest - dat_1.before_rest, dat_1.BDI_month_6_post_25mg - dat_1.BDI_Baseline);
T.diff{'Q_BDI_change', :} = [R c3nl_RCI(R,height(dat_1), 1) P];


%% Plot figure

% Set up figure elements with panel
fg1 = figure();
    p=panel();
    p.fontsize = plt.fntSze;
    p.margin = [20 15 5 5];

    p.pack('h', {0.4 0.6});
    p(2).marginleft = 40;
    p(2).pack('v', {0.45 0.5});
    p(2,1).marginbottom = 40;
    p(2,1).pack('h', {0.5 0.5});
    p(2,2).pack('h',{0.001, 0.99})

    
% Scatter histogram plot of individual's Q for each session
p(1).select(); hold on 
    offSt=0.25;
    densityScaler = 0.25;
    bwdth = 0.2;
    falpha = 0.3;

    % Plot the density for session 1
        [f,xi]=ksdensity(Q1(:,1),'Bandwidth',bwdth);
        f=-(f*densityScaler) - offSt;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'control',:});
        patch(f,xi,plt.colour{'control',:},'FaceAlpha',falpha,'Edgecolor','none');
    % plot median & mean lines
        [~,mdnI]=min(abs(median(Q1(:,1)) - xi)); plot([f(mdnI) max(f)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth)
        [~,mnI]=min(abs(mean(Q1(:,1)) - xi)); plot([f(mnI) max(f)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth)                
    % Plot the density for session 2
        [f,xi]=ksdensity(Q1(:,2),'Bandwidth',bwdth);
        f=f*densityScaler + offSt ;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'psilocybin',:});    
        patch(f,xi,plt.colour{'psilocybin',:},'FaceAlpha',falpha,'Edgecolor','none');
    % plot median & mean lines
        [~,mdnI]=min(abs(median(Q1(:,2)) - xi));plot([min(f) f(mdnI)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'psilocybin',:},'LineWidth',plt.lneWdth)
        [~,mnI]=min(abs(mean(Q1(:,2)) - xi)); plot([min(f) f(mnI)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'psilocybin',:},'LineWidth',plt.lneWdth)               
    % Plot individual data points
        [~,xyLoc]=plotSpread(Q1, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','psilocybin'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5);    
        tmpX=reshape(xyLoc(:,1), [], 2)';
    % Plot individual differences between sessions
        plot(tmpX(:,diff(Q1,[],2)>0), Q1(diff(Q1,[],2)>0,:)','Color',[0.85 0.85 0.85 0.75],'LineWidth',plt.lneWdth);
        plot(tmpX(:,diff(Q1,[],2)<0), Q1(diff(Q1,[],2)<0,:)','Color',[plt.colour{'psilocybin',:} 0.8],'LineWidth',plt.lneWdth);
        plotSpread(Q1, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','psilocybin'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5)         
    % Style
        set(gca,'Box','on','GridAlpha',0.05,'XTick',[-offSt offSt], 'XTickLabel',{'Baseline','Post-treatment'},'FontSize',plt.fntSze,'YTick', 1:4,'LineWidth', plt.lneWdth)
        ylim([0.3 4.15]); xlim([max(abs(xlim))*-1 max(abs(xlim))])
        ylabel('Brain modularity ({\itQ})');    
    
% Plot the correlation between Follow-up Q & 6 month BDI    
p(2,1,1).select(); hold on
    X=dat_1.after_rest;            
    Y=dat_1.BDI_month_6_post_25mg; 

    scatter(X,Y,plt.mrkSze*4,'k','filled'); 
    ls=lsline;ls.LineWidth=3; ls.Color=[0 0 0]; ls.LineWidth=plt.lneWdth;

    set(gca, 'FontSize', plt.fntSze*0.8, 'XTick', 1:0.5:3.5, 'LineWidth', plt.lneWdth,'YTick',0:10:55,'XTick',1:.5:4,'LineWidth',plt.lneWdth,'box','on')
    xlim([0.9 3.6]); ylim([-5 55])
    ylabel('BDI: 6 months'); xlabel('Post-treatment modularity') 

    % Add text corr stat to plot
    tmpR = num2str(round(T.corr{'after_rest:BDI_month_6_post_25mg',{'R'}},2));
    tmpP = num2str(round(T.corr{'after_rest:BDI_month_6_post_25mg',{'p_fdr'}},2));
    text(1,45,['r = ' tmpR '\newlinep = ' tmpP], 'FontSize', plt.fntSze)
    
    
% Plot correlation between difference in Q & BDI    
p(2,1,2).select(); hold on
    X=dat_1.after_rest - dat_1.before_rest;             % Neg = Q was lower after
    Y=dat_1.BDI_month_6_post_25mg - dat_1.BDI_Baseline; % Neg = BDI was lower after

    scatter(X,Y,plt.mrkSze*4,'k','filled'); 
    ls=lsline;ls.LineWidth=3; ls.Color=[0 0 0]; ls.LineWidth=plt.lneWdth;
    
    set(gca, 'YTick', min(Y):10:10, 'FontSize', plt.fntSze*0.8,'LineWidth',plt.lneWdth,'box','on')
    xlim([-1.05 0.55]); ylim([-45 15])
    ylabel('BDI change: 6 months'); xlabel('Modularity change') 
    
    % Add stat text to plot
    tmpR = num2str(round(T.diff{'Q_BDI_change',{'R'}},2));
    tmpP = num2str(round(T.diff{'Q_BDI_change',{'p'}},2));
    text(-0.98,6,['r = ' tmpR '\newlinep = ' tmpP], 'FontSize', plt.fntSze)
    
    
% Plot changes in DMN recruitment/integration     
p(2,2,2).select(); hold on
    x_lim = [0.4 3.6];
    line0 = plot([0 0], x_lim, 'Color','k','LineWidth',plt.lneWdth);
    boxplot(dat_1{:,{'DMN_R_diff','DMN_EN_I_diff','DMN_SN_I_diff'}},'Widths',0.5,'Colors',plt.colour{'control',:},'orientation','horizontal')
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1 
        h(j).LineWidth=plt.lneWdth*0.5;
        patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',1,'EdgeColor',h(j).Color);
    end 
    boxplot(dat_1{:,{'DMN_R_diff','DMN_EN_I_diff','DMN_SN_I_diff'}},'Widths',0.5,'Colors',plt.colour{'control',:},'orientation','horizontal');
    
    set(findobj(gca,'type','line'),'linew',plt.lneWdth*0.5,'Color','k');
    set(findobj(gca,'Tag','Median'),'linew',plt.lneWdth, 'Color','k')
    set(findobj(gca,'Tag','Box'),'Color','k');    
    set(gca, 'YDir','reverse','box','on', 'LineWidth',plt.lneWdth,'YTickLabel',{'DMN','DMN-EN','DMN-SN'},'YTick',-2:3)
    line0.LineWidth=plt.lneWdth;
    xlabel('Network recruitment/integration')
    
% Define final figure size
set(fg1, 'Position', fgPos)

% Save out SVG & jpeg
if saveFigs
    print('../figures/psilodep1_Q_figure.svg', '-dsvg')
    saveas(fg1, '../figures/psilodep1_Q_figure.jpg')
end



%% PSILODEP 2 %%

fdrSz = 50; % FDR-corrected star '*' text font size
lwdth=3; % Width of plot Box line
lbl = flex.R.psilocybin.Row; 
rTicks = [0:-0.2:-1]'; % Corr axis label

offSt=0.25;
densityScaler = 0.25;
bwdth = 0.2;
falpha = 0.3;

fg2 = figure();
    p=panel();
    p.fontsize = plt.fntSze;
    p.margin = [20 25 20 10];

    % Pack fig spaces
    p.pack('h', {0.5,[],[],0.03});
    p(1).pack('v', {0.48 0.48}); 
    p(2).pack('v', {0.48 0.48}); 
    p(3).pack('v', {0.48 0.48}); 
    p(4).pack('v', {0.48 0.48}); 
    
    % Add in margin
     p(1).marginright = 45;
     p(2).marginright = 30;
     p(3).marginright = 0;
     p(1,1).marginbottom = 50;
     p(2,1).marginbottom = 50;
     p(3,1).marginbottom = 50;
     p(4,1).marginbottom = 50;

p(1,1).select();
% psilocybin
    tmp = dat_2{dat_2.arm=='psilocybin',{'ses_1','ses_2'}};
    hold on;
    % Plot the density for session 1
        [f,xi]=ksdensity(tmp(:,1),'Bandwidth',bwdth);
        f=-(f*densityScaler) - offSt;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'control',:});
        patch(f,xi,plt.colour{'control',:},'FaceAlpha',falpha,'Edgecolor','none');
        % plot median & mean lines
        [~,mdnI]=min(abs(median(tmp(:,1)) - xi)); plot([f(mdnI) max(f)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth*0.5)
        [~,mnI]=min(abs(mean(tmp(:,1)) - xi)); plot([f(mnI) max(f)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth*0.5)                   
    % Plot the density for session 2
        [f,xi]=ksdensity(tmp(:,2),'Bandwidth',bwdth);
        f=f*densityScaler + offSt ;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'psilocybin',:});    
        patch(f,xi,plt.colour{'psilocybin',:},'FaceAlpha',falpha,'Edgecolor','none');
        % plot median & mean lines
        [~,mdnI]=min(abs(median(tmp(:,2)) - xi));plot([min(f) f(mdnI)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'psilocybin',:},'LineWidth',plt.lneWdth*0.5)
        [~,mnI]=min(abs(mean(tmp(:,2)) - xi)); plot([min(f) f(mnI)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'psilocybin',:},'LineWidth',plt.lneWdth*0.5)             
    % Plot individual data points
        [~,xyLoc]=plotSpread(tmp, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','psilocybin'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5);    
        tmpX=reshape(xyLoc(:,1), [], 2)';
    % Plot individual differences between sessions
        plot(tmpX(:,diff(tmp,[],2)>0), tmp(diff(tmp,[],2)>0,:)','Color',[0.85 0.85 0.85 0.75],'LineWidth',plt.lneWdth*0.5);
        plot(tmpX(:,diff(tmp,[],2)<0), tmp(diff(tmp,[],2)<0,:)','Color',[plt.colour{'psilocybin',:} 0.8],'LineWidth',plt.lneWdth*0.75);
        plotSpread(tmp, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','psilocybin'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5)         
    % Style
        set(gca,'Box','on','GridAlpha',0.05,'XTick',[-offSt offSt], 'XTickLabel',{'Baseline','Post-treatment'},'FontSize',plt.fntSze,'YTick', 1:4,'LineWidth', plt.lneWdth)
        ylim([0.2 4.25]); xlim([max(abs(xlim))*-1 max(abs(xlim))])
        ylabel('Brain modularity ({\itQ})'); 
    title('Psilocybin') 


p(1,2).select();
    % escitalopram
    tmp = dat_2{dat_2.arm=='escitalopram',{'ses_1','ses_2'}};
    hold on
    % Plot the density for session 1
        [f,xi]=ksdensity(tmp(:,1),'Bandwidth',bwdth);
        f=-(f*densityScaler) - offSt;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'control',:});
        patch(f,xi,plt.colour{'control',:},'FaceAlpha',falpha,'Edgecolor','none');
        % plot median & mean lines
        [~,mdnI]=min(abs(median(tmp(:,1)) - xi)); plot([f(mdnI) max(f)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth*0.5)
        [~,mnI]=min(abs(mean(tmp(:,1)) - xi)); plot([f(mnI) max(f)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'control',:},'LineWidth',plt.lneWdth*0.5)                
    % Plot the density for session 2
        [f,xi]=ksdensity(tmp(:,2),'Bandwidth',bwdth);
        f=f*densityScaler + offSt ;
        plot(f,xi,'LineWidth',plt.lneWdth,'Color',plt.colour{'escitalopram',:});    
        patch(f,xi,plt.colour{'escitalopram',:},'FaceAlpha',falpha,'Edgecolor','none');
        % plot median & mean lines
        [~,mdnI]=min(abs(median(tmp(:,2)) - xi));plot([min(f) f(mdnI)],repmat(xi(mdnI),1,2),':','MarkerSize',plt.mrkSze,'Color',plt.colour{'escitalopram',:},'LineWidth',plt.lneWdth*0.5)
        [~,mnI]=min(abs(mean(tmp(:,2)) - xi)); plot([min(f) f(mnI)],repmat(xi(mnI),1,2),'-','MarkerSize',plt.mrkSze,'Color',plt.colour{'escitalopram',:},'LineWidth',plt.lneWdth*0.5)               
    % Plot individual data points
        [~,xyLoc]=plotSpread(tmp, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','escitalopram'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5);    
        tmpX=reshape(xyLoc(:,1), [], 2)';
    % Plot individual differences between sessions
        plot(tmpX(:,diff(tmp,[],2)>0), tmp(diff(tmp,[],2)>0,:)','Color',[0.85 0.85 0.85 0.75],'LineWidth',plt.lneWdth*0.5);
        plot(tmpX(:,diff(tmp,[],2)<0), tmp(diff(tmp,[],2)<0,:)','Color',[plt.colour{'escitalopram',:} 0.8],'LineWidth',plt.lneWdth*0.75);
        plotSpread(tmp, 'xValues',[-offSt offSt],'distributionColors',plt.colour{{'control','escitalopram'},:},'MarkerSize',plt.mrkSze, 'spreadWidth',0.5)         
    % Style
        set(gca,'Box','on','GridAlpha',0.05,'XTick',[-offSt offSt], 'XTickLabel',{'Baseline','Post-treatment'},'FontSize',plt.fntSze,'YTick', 1:4,'LineWidth', plt.lneWdth)
        ylim([0.2 4.25]); xlim([max(abs(xlim))*-1 max(abs(xlim))])
        ylabel('Brain modularity ({\itQ})');  
    title('Escitalopram')       

    
    
% plot Q & BDI 6 week changes
p(2,1).select();
    X=diff(dat_2{dat_2.arm=='psilocybin', {'ses_1','ses_2'}},[],2); % Neg = reduced Q
    Y=diff(dat_2{dat_2.arm=='psilocybin', {'BDI_BASELINE1wpredd1','BDI_6WEEKSpostdd1'}},[],2); % Neg = reduced depression
    h=scatter(X,Y,100,'k','filled'); hold on
    ls=lsline;ls.LineWidth=3; ls.Color=[0 0 0]; ls.LineWidth=plt.lneWdth;
    set(gca, 'FontSize', plt.fntSze,'LineWidth',plt.lneWdth,'box','on')%,'YTick',0:5:55,'XTick',1:.5:4)    
    ylabel('BDI change: 6 weeks', 'FontSize', plt.fntSze)
    xlabel('Modularity change', 'FontSize', plt.fntSze) 
    ylim([-42 12]); xlim([-2.4 1.4]);
    
    % Add text corr stat to plot
    [R,P]=corr(X,Y, 'tail','right'); %1 tailed as informed by study1
    tmpR = num2str(round(R,2)); tmpP = num2str(round(P,3));
    text(-2.2,9.5,['r = ' tmpR], 'FontSize', plt.fntSze*0.8)
    text(-2.2,4.5,['p = ' tmpP], 'FontSize', plt.fntSze*0.8)   

p(2,2).select();
    X=diff(dat_2{dat_2.arm=='escitalopram', {'ses_1','ses_2'}},[],2); % Neg = reduced Q
    Y=diff(dat_2{dat_2.arm=='escitalopram', {'BDI_BASELINE1wpredd1','BDI_6WEEKSpostdd1'}},[],2); % Neg = reduced depression
    h=scatter(X,Y,100,'k','filled'); hold on
    ls=lsline;ls.LineWidth=3; ls.Color=[0 0 0]; ls.LineWidth=plt.lneWdth;
    set(gca, 'FontSize', plt.fntSze,'LineWidth',plt.lneWdth,'box','on')%,'YTick',0:5:55,'XTick',1:.5:4)    
    ylabel('BDI change: 6 weeks', 'FontSize', plt.fntSze)
    xlabel('Modularity change', 'FontSize', plt.fntSze) 
    ylim([-42 12]); xlim([-2.4 2.4]);
    
    % Add text corr stat to plot
    [R,P]=corr(X,Y, 'tail','right'); %1 tailed as informed by study1
    tmpR = num2str(round(R,2)); tmpP = num2str(round(P,3));
    text(-2.2,9.5,['r = ' tmpR], 'FontSize', plt.fntSze*0.8)
    text(-2.2,4.5,['p = ' tmpP], 'FontSize', plt.fntSze*0.8)  
    
    
a=axes;      
    h=imagesc(flex.R.psilocybin{:,:}, [-1 1]); set(h, 'AlphaData', flex.P.psilocybin{:,:}<=0.05); % Plot R values and threshold visualisation with unc p<0.05
    axis('tight')
    set(gca,'box','on','YDir','reverse','YTickLabel', lbl, 'YTick', 1:numel(lbl), 'XTick',1:numel(lbl), 'XTickLabel', lbl, 'LineWidth', plt.lneWdth,'FontSize',plt.fntSze); 
    colormap(a, fWS_cmap([plt.colour{'psilocybin',:};0.95 0.95 0.95]))
    
    % Add significance stars to those corrs that survive FDR-correction
    tmpIx = tril(true(numel(lbl)));
    tmpP = table2array(flex.P.psilocybin);
    tmp=zeros(numel(lbl));
    [~, ~, ~, Pc]=fdr_bh(tmpP(tmpIx), 0.05, 'pdep', 'no');
    tmp(tmpIx) = Pc<=0.05; 
    [tmpI,tmpJ]=ind2sub(size(tmp), find(tmp)); 
    text([tmpJ; tmpI], [tmpI; tmpJ]+0.25,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fdrSz); 
p(3,1).select(a);
    

a=axes;      
    h=imagesc(flex.R.escitalopram{:,:}, [-1 1]); set(h, 'AlphaData', flex.P.escitalopram{:,:}<=0.05); % Plot R values and threshold visualisation with unc p<0.05
    axis('tight')
    set(gca,'box','on','YDir','reverse','YTickLabel', lbl, 'YTick', 1:numel(lbl), 'XTick',1:numel(lbl), 'XTickLabel', lbl, 'LineWidth', plt.lneWdth,'FontSize',plt.fntSze); 
    colormap(a, fWS_cmap([plt.colour{'escitalopram',:};0.95 0.95 0.95]))
    
    % Add significance stars to those corrs that survive FDR-correction
    tmpIx = tril(true(numel(lbl)));
    tmpP = table2array(flex.P.escitalopram);
    tmp=zeros(numel(lbl));
    [~, ~, ~, Pc]=fdr_bh(tmpP(tmpIx), 0.05, 'pdep', 'no');
    tmp(tmpIx) = Pc<=0.05; 
    [tmpI,tmpJ]=ind2sub(size(tmp), find(tmp)); 
    text([tmpJ; tmpI], [tmpI; tmpJ],'*','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fdrSz); 
p(3,2).select(a);        
    
    
a=axes; % plot corr colour bar
    imagesc(a, rTicks(rTicks<=0)); 
    set(gca, 'TickLength',[0 0],'YAxisLocation','right','YTickLabel',strrep(strsplit(num2str(rTicks')), '0.','.'),'YTick',1:numel(rTicks),'XTick',[],'LineWidth',lwdth*0.75);
    axis('tight')
    title(a,'{\it r}','FontSize',plt.fntSze)    
    colormap(a, fWS_cmap([plt.colour{'psilocybin',:};0.95 0.95 0.95]))
p(4,1).select(a);
    

set(fg2, 'Position', fgPos)

% Save out SVG & jpeg
if saveFigs
    print('../figures/psilodep2_Q_figure.svg', '-dsvg')
    saveas(fg2, '../figures/psilodep2_Q_figure.jpg')
end
