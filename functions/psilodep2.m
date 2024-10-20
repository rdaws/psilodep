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
%
% Psilodep 2:- Preprocessed data handling & visulisation
% 
% Requires 
%    panel.m (https://uk.mathworks.com/matlabcentral/fileexchange/20003-panel)
%
%% Load in data
clear; close all

% load in data from each study & plotting variables
load('../data/psilodep1/dat_1.mat');
load('../data/psilodep2/dat_2.mat');
load('../data/psilodep2/flex.mat');
load('../plotting_vars.mat');
    

%% Between trial baseline BDI

% Stat table (T) column labels & row name 
    stLbl = {'mean_diff','t','CI','p','cohenD'}; 
    cpLbl = {'baseline_BDI_p1p2'};
    
    T.baseline_BDI = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.baseline_BDI.CI = NaN(height(T.baseline_BDI), 2);
    [~, P, CI, tST] = ttest2(dat_1.BDI_Baseline, dat_2.BDI_BASELINE1wpredd1);
    
    tmp = [dat_1.BDI_Baseline; dat_2.BDI_BASELINE1wpredd1];
   
    T.baseline_BDI{:,:} = [mean(dat_1.BDI_Baseline)-mean(dat_2.BDI_BASELINE1wpredd1) tST.tstat CI' P c3nl_cohens_D(tmp, [ones(height(dat_1),1); 2*ones(height(dat_2),1)])];
    
    

%% BDI changes & between arm comparison
    
% Extract BDI variables names
    tmplbl=dat_2.Properties.VariableNames(contains(dat_2.Properties.VariableNames, 'BDI'));
    M=[char(join(tmplbl,',')) ' ~ arm + 1']; % RM with treament arm as a between factor
    
    W=table(categorical((1:numel(tmplbl))'), 'VariableNames', {'timepoint'});
    AT=ranova(fitrm(dat_2, M, 'WithinDesign', W),'WithinModel','timepoint');

% Stat table (T) column labels 
    stLbl = {'mean_diff','t','CI','p','cohenD'}; 
% Comparison labels of BDI timepoints seperated by a ":"
    cpLbl = {'BDI_BASELINE1wpredd1:BDI_2weekspostdd1', 'BDI_BASELINE1wpredd1:BDI_4weekspostdd1', 'BDI_BASELINE1wpredd1:BDI_6WEEKSpostdd1'};

% Create stat table
    T.BDI = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.BDI.CI = NaN(height(T.BDI), 2);

% For each comparison, calc t-stats & effect size
    for ii = 1:numel(cpLbl)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.BDI.Row(ii),':'); 
        % Two sample ttest on the baseline diff
        tmpP = diff(dat_2{dat_2.arm=='psilocybin', tmp},[],2);
        tmpE = diff(dat_2{dat_2.arm=='escitalopram', tmp},[],2);
        
        [~,P,CI,tST]=ttest2(tmpP,tmpE);
        % put tstat mean diff and cohen's D into the table 
        T.BDI{T.BDI.Row(ii), :} = [mean(tmpP)-mean(tmpE) tST.tstat CI' P c3nl_cohens_D([tmpP;tmpE], [ones(numel(tmpP),1); 2*ones(numel(tmpE),1)])];
    end    
    [~, ~, ~, Pc]=fdr_bh(T.BDI.p, 0.05, 'pdep', 'no');
    T.BDI.p_fdr = Pc;

    
%% Modularity session difference for each arm

% Stat table (T) column labels 
    stLbl = {'mean_diff','t','CI','p','cohenD'}; 
% Comparison labels of BDI timepoints seperated by a ":"
    cpLbl = {'psilocybin:ses_1:ses_2', 'escitalopram:ses_1:ses_2'};

% Create stat table
    T.Q = array2table(NaN(numel(cpLbl),numel(stLbl)), 'VariableNames', stLbl,'RowNames',cpLbl);
    T.Q.CI = NaN(height(T.Q), 2);

% For each comparison, calc t-stats & effect size
    for ii = 1:numel(cpLbl)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.Q.Row(ii),':'); 
        % Two sample ttest on the baseline diff
        tmpQ = diff(dat_2{dat_2.arm==tmp{1}, tmp(2:3)},[],2);
        
        [~,P,CI,tST]=ttest(tmpQ);
        % put tstat mean diff and cohen's D into the table 
        T.Q{T.Q.Row(ii), :} = [mean(tmpQ) tST.tstat CI' P c3nl_cohens_D(tmpQ)];
    end        


%% Modularity change stats and BDI correlations     
       
    tmplbl=dat_2.Properties.VariableNames(contains(dat_2.Properties.VariableNames, 'ses_'));
    M=[char(join(tmplbl,',')) ' ~ arm + 1'];
    W=table(categorical((1:numel(tmplbl))'), 'VariableNames', {'timepoint'});
    AT=ranova(fitrm(dat_2, M, 'WithinDesign', W),'WithinModel','timepoint');
    
% Preallocate modularity table    
    cpLbl = {'ses_1:ses_2'}; 
    grp = {'psilocybin','escitalopram'};
    stLbl = {'mean_diff','t','CI','p','cohenD'};    
    T.Q = array2table(NaN(numel(cpLbl)*numel(grp),numel(stLbl)), 'VariableNames', stLbl,'RowNames',strcat(grp', ':', cpLbl));
    T.Q.CI = NaN(height(T.Q), 2);
    
% Effect of psilocybin-therapy on resting-state Q
    for ii = 1:height(T.Q)
        for jj = 1:numel(grp)
            % Split the string with ":" delim to indentify comparison timepoints
            tmp = split(T.Q.Row(ii),':'); 
            % Paired ttest
            [~,P,CI,tST]=ttest(dat_2{dat_2.arm==tmp{1}, tmp(2)}, dat_2{dat_2.arm==tmp{1}, tmp(3)});
            % put tstat mean diff and cohen's D into the table 
            T.Q{T.Q.Row(ii), :} = [mean(diff(dat_2{dat_2.arm==tmp{1},tmp(2:3)},[],2)) tST.tstat CI' P c3nl_cohens_D(dat_2{dat_2.arm==tmp{1},tmp(2:3)})];
        end
    end
     
    
    
%% Correlation between the changes in BDI and Q (expected low Q relating to low BDI)  
    stLbl = {'R','CI','p'}; 
    grp = {'psilocybin','escitalopram'};
    cpLbl = {'BDI_BASELINE1wpredd1:BDI_6WEEKSpostdd1'};
    T.diff = array2table(NaN(numel(cpLbl)*numel(grp),numel(stLbl)), 'VariableNames', stLbl,'RowNames',strcat(repelem(grp',numel(cpLbl),1), ':', repmat(cpLbl',numel(grp),1)));
    T.diff.CI=NaN(height(T.diff), 2);
    
% Relationship between 
    for ii = 1:height(T.diff)
        % Split the string with ":" delim to indentify comparison timepoints
        tmp = split(T.diff.Row(ii),':'); 
        % Pearson corr - 1- tailed
        [R,P]=corr(diff(dat_2{dat_2.arm==tmp{1}, {'ses_1','ses_2'}},[],2), diff(dat_2{dat_2.arm==tmp{1}, tmp(2:3)},[],2),'tail','right'); 
        % put effects in table
        T.diff{T.diff.Row(ii), :} = [R c3nl_RCI(R,sum(dat_2.arm==tmp{1})) P];
    end  
