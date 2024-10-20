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
% Psilodep :- Stats reported in manuscript
% 
% Requires 
%    third_party folder and subfolders to be in path 

%%
clear; close all

% load in data & plotting variables
load('../data/psilodep1/dat_1.mat');
load('../plotting_vars.mat');


%% Beck dession inventory stats

% Stat table (T) column labels 
    stLbl = {'mean_diff','t','CI','p','cohenD'}; 
% Comparison labels of BDI timepoints seperated by a ":"
    cpLbl = {'BDI_Baseline:BDI_week_1_post_25mg','BDI_Baseline:BDI_month_6_post_25mg'};

% Create stat table
    T.BDI = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.BDI.CI = NaN(height(T.BDI), 2);

% For each comparison, calc t-stats & effect size
    for ii = 1:numel(cpLbl)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.BDI.Row(ii),':'); 
        % Paired ttest
        [~,P,CI,tST]=ttest(dat_1{:, tmp(1)}, dat_1{:, tmp(2)});
        % put tstat mean diff and cohen's D into the table 
        T.BDI{T.BDI.Row(ii), :} = [mean(diff(dat_1{:,tmp},[],2)) tST.tstat CI' P c3nl_cohens_D(dat_1{:,tmp})];
    end
    
    
%% modularity change

%
% Pearson FC from 100 ROIs (Schaefer et al., 2018), fisher tranformed and positive values retained.
% Modularity genLouvain.m estimated 100 times and partition with largest Q normalised against the mean Q from 100 runs of randomly rewired (shuffled) FC matrices
%

% Preallocate modularity table    
    cpLbl = {'before_rest:after_rest'}; 
    T.Q = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.Q.CI = NaN(height(T.Q), 2);

% Effect of psilocybin-therapy on resting-state Q
    for ii = 1:numel(cpLbl)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.Q.Row(ii),':'); 
        % Paired ttest
        [~,P,CI,tST]=ttest(dat_1{:, tmp(1)}, dat_1{:, tmp(2)});
        % put tstat mean diff and cohen's D into the table 
        T.Q{T.Q.Row(ii), :} = [mean(diff(dat_1{:,tmp},[],2)) tST.tstat CI' P c3nl_cohens_D(dat_1{:,tmp})];
    end

    
%% BDI modularity correlations (expected low Q relating to low BDI)
    cpLbl = {'after_rest:BDI_week_1_post_25mg','after_rest:BDI_month_3_post_25mg','after_rest:BDI_month_6_post_25mg'}; 
    stLbl = {'R','CI','p'};
    T.corr = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.corr.CI = NaN(height(T.corr), 2);

    for ii = 1:numel(cpLbl)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.corr.Row(ii),':'); 
        % One-tailed pearson correlation
        [R,P] = corr(dat_1{:, tmp(1)}, dat_1{:, tmp(2)});
        % put tstat mean diff and cohen's D into the table 
        T.corr{T.corr.Row(ii), :} = [R c3nl_RCI(R,height(dat_1),1) P];
    end    
    
% Add in FDR correction
    [~, ~, ~, Pc]=fdr_bh(T.corr.p, 0.05, 'pdep', 'no');
    T.corr.p_fdr = Pc;
    
% Correlation between the changes in BDI and Q (expected low Q relating to low BDI)  
    T.diff = array2table(NaN(1,numel(stLbl)), 'VariableNames', stLbl,'RowNames',{'Q_BDI_change'});
    T.diff.CI=NaN(height(T.diff), 2);
    
    [R,P] = corr(dat_1.after_rest - dat_1.before_rest, dat_1.BDI_month_6_post_25mg - dat_1.BDI_Baseline);
    T.diff{'Q_BDI_change', :} = [R c3nl_RCI(R,height(dat_1), 1) P];
    
    
%% Recruitment / integration analysis       

% Stat table (T) column labels 
    stLbl = {'mean_diff','t','CI','p','cohenD'}; 
    cpLbl = {'DMN','DMN_EN','DMN_SN'};
% Create stat table
    [T.I, T.R] = deal(array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl));
    [T.I.CI, T.R.CI] = deal(NaN(height(T.R), 2));
    
% A priori depression network tests    
    % DMN recruitment down
    [~,tmpP,CI,tmpT] = ttest(dat_1.DMN_R_diff);
    T.R{'DMN',:} = [mean(dat_1.DMN_R_diff), tmpT.tstat CI' tmpP c3nl_cohens_D(dat_1.DMN_R_diff)];
    
    % DMN - FP integration up
    [~,tmpP,CI,tmpT] = ttest(dat_1.DMN_EN_I_diff);
    T.I{'DMN_EN',:} = [mean(dat_1.DMN_EN_I_diff), tmpT.tstat CI' tmpP c3nl_cohens_D(dat_1.DMN_EN_I_diff)];
    
    % DMN - FP integration up
    [~,tmpP,CI,tmpT] = ttest(dat_1.DMN_SN_I_diff);
    T.I{'DMN_SN',:} = [mean(dat_1.DMN_SN_I_diff), tmpT.tstat CI' tmpP c3nl_cohens_D(dat_1.DMN_SN_I_diff)];
    % Add in fdr p
    [~, ~, ~, Pc]=fdr_bh(T.I.p(2:3), 0.05, 'pdep', 'no');
    T.I.p_fdr = [NaN; Pc];    
        