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
% Psilodep :- Beck's Depression Inventory (BDI) visualisation
% 
% Requires 
%    third_party folder and subfolders to be in path 

%%
clear; close all

saveFigs = true;        % Boolean to save out figs or not.
fgPos = [1 1 898 796];  % Final figure size

% load in data from each study & plotting variables
load('../data/psilodep1/dat_1.mat');
load('../data/psilodep2/dat_2.mat');
load('../plotting_vars.mat');

clim = [14 28];     % colour map range for individual BDI raster plots
a = [0.5, 1];       % Alpha of boxes for baseline study1 vs study2 BDI

% Extract BDI measures for each study
BDI_1 = dat_1{:,contains(dat_1.Properties.VariableNames, 'BDI')};
xlb_1 = {'BL','1w','3m','6m'};

BDI_2 = dat_2{:,contains(dat_2.Properties.VariableNames, 'BDI')};
xlb_2 = {'BL','2w','4w','6w'};

%% PLOTTING

% Create figure elements with panel
fg = figure();
    p=panel();
    p.fontsize = plt.fntSze;
    p.margin = [20 20 5 5];

    % Pack fig spaces
    p.pack('h', {0.7, 0.23, 0.015}); % Major columns
    p(1).pack('v', {0.5, 0.5});      % Space for the BDI timecourse boxplots 
    
    % Spaces for individual raster plots
    p(2).pack('v', {0.5, 0.5}); p(2).marginleft = 20; 
    p(2,1).pack('v', {0.5, 0.5}); 
    p(2,2).pack('v', {0.5, 0.5}); 
    
    % Spaces for colour bars
    p(3).pack('v', {0.5, 0.5}); p(3).marginleft=5;
    p(3,1).pack('v', {0.5, 0.5})

% Psilodep1 BDI timecourse boxplot
p(1,1).select(); hold on
    boxplot(BDI_1, 'Widths', 0.4, 'Colors', plt.colour{'psilocybin',:});
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1 
        h(j).LineWidth=plt.lneWdth*0.5;
        patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',0.8,'EdgeColor',h(j).Color);
    end 
    boxplot(BDI_1, 'Widths', 0.4, 'Colors', plt.colour{'psilocybin',:});
    
    ylabel('BDI')
    set(findobj(gca,'type','line'),'linew',plt.lneWdth*0.5,'Color','k');
    set(findobj(gca,'Tag','Median'),'linew',plt.lneWdth, 'Color','k'); set(findobj(gca,'Tag','Box'),'Color','k'); 
    set(gca,'XTickLabel', xlb_1, 'XGrid','On','YGrid','On', 'GridColor',[1 1 1], 'Box','on','LineWidth',plt.lneWdth, 'FontSize',plt.fntSze*1.5)    
    
% Psilodep1 individual BDI raster    
imP = p(2,1,1).select();
    [~,I]=sort(sum(BDI_1,2));
    imagesc(flipud(BDI_1(I,:)),clim);  
    set(gca, 'XTickLabel', xlb_1, 'XTick', 1:numel(xlb_1), 'YTick', [], 'FontSize', plt.fntSze, 'Box','On','LineWidth',plt.lneWdth)
    colormap(imP, 'gray'); colormap(imP, flipud(imP.Colormap));imP.ColorScale = 'linear';
    xlim([0.5 size(BDI_1,2) + 0.5]); ylim([0.5 size(BDI_1,1) + 0.5])    

% Psilodep1 individual BDI colour bar
cbP = p(3,1,1).select();    
    imagesc(linspace(clim(2), clim(1))')  
    set(gca,'Box','On', 'YAxisLocation', 'right', 'XTick',[], 'YTick', [15 75], 'YTickLabel', {'Mild','Severe'},'YTickLabelRotation',90,'FontSize',plt.fntSze*0.75,'LineWidth',plt.lneWdth*0.5,'TickLength',[0 0])
    axis('tight'); colormap(cbP, 'gray');    

% Psilodep1/2 baseline BDI comparison    
p(2,1,2).select(); hold on
    boxplot([dat_1.BDI_Baseline; dat_2.BDI_BASELINE1wpredd1], [zeros(height(dat_1),1); ones(height(dat_2),1)],'Widths',0.3,'Colors',plt.colour{'control',:});
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1 
        h(j).LineWidth=plt.lneWdth*0.5;
        patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',a(j),'EdgeColor',h(j).Color);
    end 
    boxplot([dat_1.BDI_Baseline; dat_2.BDI_BASELINE1wpredd1], [zeros(height(dat_1),1); ones(height(dat_2),1)],'Widths',0.3,'Colors',plt.colour{'control',:});   
    
    ylabel('Baseline BDI','FontSize',plt.fntSze)
    set(findobj(gca,'type','line'),'linew',plt.lneWdth*0.5,'Color','k');
    set(findobj(gca,'Tag','Median'),'linew',plt.lneWdth, 'Color','k')
    set(findobj(gca,'Tag','Box'),'Color','k'); 

    set(gca,'XTickLabel',{'TRD','MDD'},'YTick',[0 25 50],'YAxisLocation','right', 'GridColor',[1 1 1], 'Box','on','LineWidth',plt.lneWdth, 'FontSize',plt.fntSze)    
    ylim(p(1,1).axis.YLim)    
    
% Psilodep2 BDI timecourse boxplot    
p(1,2).select(); hold on
    xPos = [1:2 4:5 7:8 10:11];
    boxplot(BDI_2,[repelem((1:size(BDI_2,2))',size(BDI_2,1),1) double(repmat(dat_2.arm,size(BDI_2,2),1))],...
            'Widths',0.8,'Colors',repmat(plt.colour{{'escitalopram','psilocybin'},:},size(BDI_2,2),1),'Positions',xPos);
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1 
        h(j).LineWidth=plt.lneWdth*0.5;
        patch(get(h(j),'Xdata'),get(h(j),'Ydata'),h(j).Color,'FaceAlpha',0.8,'EdgeColor',h(j).Color);
    end 
    boxplot(BDI_2,[repelem((1:size(BDI_2,2))',size(BDI_2,1),1) double(repmat(dat_2.arm,size(BDI_2,2),1))],...
            'Widths',0.8,'Colors',repmat(plt.colour{{'escitalopram','psilocybin'},:},size(BDI_2,2),1),'Positions',xPos);
       
    ylabel('BDI','FontSize',plt.fntSze)
    set(findobj(gca,'type','line'),'linew',plt.lneWdth*0.5,'Color','k');
    set(findobj(gca,'Tag','Median'),'linew',plt.lneWdth, 'Color','k')
    set(findobj(gca,'Tag','Box'),'Color','k'); 
    ylim(p(1,1).axis.YLim); xlim([0 12])
    
    set(gca, 'XTick',xPos,'XGrid','On','YGrid','On', 'GridColor',[1 1 1], 'Box','on','LineWidth',plt.lneWdth, 'FontSize',plt.fntSze)    
    legend({'Escitalopram','Psilocybin'},'Box','Off','FontSize',plt.fntSze);    
    
    % Manual xtick labels
    xPos = reshape(xPos,2,[])';
    for ii = 1:numel(xlb_2)
        text(mean(xPos(ii,:)),  -5, xlb_2{ii},'HorizontalAlignment','center','FontSize',plt.fntSze)
    end    
    
% Psilodep2 psilocybin arm individual raster
imP = p(2,2,1).select();
    tmpP = BDI_2(dat_2.arm=='psilocybin',:);
    [~,I] = sort(nansum(tmpP,2));
    imagesc(flipud(tmpP(I,:)),clim);
    set(gca, 'XTickLabel', xlb_2, 'XTick', 1:numel(xlb_2), 'YTick', [], 'FontSize', plt.fntSze, 'Box','on','LineWidth',plt.lneWdth)
    xlim([0.5 size(tmpP,2) + 0.5]); ylim([0.5 size(tmpP,1) + 0.5])     
% Psilodep2 escitalopram arm individual raster    
imE = p(2,2,2).select();
    tmpE = BDI_2(dat_2.arm=='escitalopram',:);
    [~,I] = sort(nansum(tmpE,2));
    imagesc(flipud(tmpE(I,:)),clim);
    set(gca, 'XTickLabel', xlb_2, 'XTick', 1:numel(xlb_2), 'YTick', [], 'FontSize', plt.fntSze, 'Box','on','LineWidth',plt.lneWdth)
    xlim([0.5 size(tmpE,2) + 0.5]); ylim([0.5 size(tmpE,1) + 0.5])             
% Psilodep2 colour bar
cbP = p(3,2).select();    
    imagesc(linspace(clim(2), clim(1))') 
    set(gca,'Box','On', 'YAxisLocation', 'right', 'XTick',[], 'YTick', [7 90], 'YTickLabel', {'Mild','Severe'},'YTickLabelRotation',90,'FontSize',plt.fntSze*0.75,'LineWidth',plt.lneWdth*0.5,'TickLength',[0 0])
    axis('tight'); colormap(cbP, 'gray');  

colormap(imP, 'gray'); colormap(imP, flipud(imP.Colormap)); imP.ColorScale = 'linear';
colormap(imE, 'gray'); colormap(imE, flipud(imE.Colormap)); imE.ColorScale = 'linear';
   
% Define final figure size
set(fg, 'Position', fgPos)

% Save out SVG & jpeg
if saveFigs
    print([Rt '/figures/BDI_figure.svg'], '-dsvg')
    saveas(fg, [Rt '/figures/BDI_figure.jpg'])
end




%% Research briefing fig


fg=figure();
    h=shadedErrorBar(1:4,mean(dat_1{:,contains(dat_1.Properties.VariableNames,'BDI')}), std(dat_1{:,contains(dat_1.Properties.VariableNames,'BDI')})/sqrt(height(dat_1)));
    h.mainLine.Color = plt.colour{'psilocybin',:};
    h.mainLine.LineWidth = plt.lneWdth*3;
    h.patch.FaceColor = plt.colour{'psilocybin',:};
    h.patch.FaceAlpha = 0.4;
    h.edge(1).Color = [plt.colour{'psilocybin',:} 0];
    h.edge(2).Color = [plt.colour{'psilocybin',:} 0];
    
    xlim([0.9 4.1]); ylim([-1 41])
    ylabel('Depression severity')
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Baseline', '1 week', '3 months', '6 months'}, 'FontSize', plt.fntSze,'LineWidth',plt.lneWdth)
   
if saveFigs
    print([Rt '/figures/BDI_errorbar_1.svg'], '-dsvg')
    saveas(fg, [Rt '/figures/BDI_errorbar_1.jpg'])
end

fg=figure();
    h=shadedErrorBar(1:4,mean(dat_2{dat_2.arm=='psilocybin',contains(dat_2.Properties.VariableNames,'BDI')}), std(dat_2{dat_2.arm=='psilocybin',contains(dat_2.Properties.VariableNames,'BDI')})/sqrt(sum(dat_2.arm=='psilocybin')));
    h.mainLine.Color = plt.colour{'psilocybin',:};
    h.mainLine.LineWidth = plt.lneWdth*3;
    h.patch.FaceColor = plt.colour{'psilocybin',:};
    h.patch.FaceAlpha = 0.4;
    h.edge(1).Color = [plt.colour{'psilocybin',:} 0];
    h.edge(2).Color = [plt.colour{'psilocybin',:} 0];
    
    
    h=shadedErrorBar(1:4,mean(dat_2{dat_2.arm=='escitalopram',contains(dat_2.Properties.VariableNames,'BDI')}), std(dat_2{dat_2.arm=='escitalopram',contains(dat_2.Properties.VariableNames,'BDI')})/sqrt(sum(dat_2.arm=='escitalopram')));
    h.mainLine.Color = plt.colour{'escitalopram',:};
    h.mainLine.LineWidth = plt.lneWdth*3;
    h.patch.FaceColor = plt.colour{'escitalopram',:};
    h.patch.FaceAlpha = 0.4;
    h.edge(1).Color = [plt.colour{'escitalopram',:} 0];
    h.edge(2).Color = [plt.colour{'escitalopram',:} 0];
    
    xlim([0.9 4.1]); ylim([-1 41])
    ylabel('Depression severity')
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Baseline', '2 weeks','4 weeks','6 weeks'}, 'FontSize', plt.fntSze,'LineWidth',plt.lneWdth)

if saveFigs
    print('../figures/BDI_errorbar_2.svg', '-dsvg')
    saveas(fg, '../figures/BDI_errorbar_2.jpg')
end    