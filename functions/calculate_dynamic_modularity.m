function [Q, modules] = calculate_dynamic_modularity(imaging_volume, atlas_volume)
%                                 _---~~(~~-_.
%                               _{        )   )
% ██████  ███████ ██████      ,   ) -~~- ( ,-' )_
% ██   ██ ██      ██   ██    (  `-,_..`., )-- '_,)
% ██████  █████   ██   ██   ( ` _)  (  -~( -_ `,  }
% ██   ██ ██      ██   ██   (_-  _  ~_-~~~~`,  ,' )
% ██   ██ ███████ ██████      `~ -^(    __;-,((()))
% Richard E. Daws  JUN2020           ~~~~ {_ -_(())
%                                          `\  }
%                                            { }   
% Psilodep :- Calculate dynamic modularity (Q)
% 
% Requires 
%    third_party folder and subfolders to be in path 
%
%
%   Args:
%       imagine_volume - 4D matrix of fMRI imaging data
%       atlas_volume   - 3D matrix of unique values denoting Regions of
%                        Interest(ROI)
% 
% Adapted from Karolina Finc et al 2020, Nat Comms. 
%
%

%% Parameters

% window size
W=15; 
gamma = 1;
omega = 1;
nRep = 100;


%% 1 Wrangle atlas

    % Find number of ROIs
    nRoi=numel(unique(atlas_volume(atlas_volume>0)));
    % Reshape atlas to 2D
    atlas_volume=reshape(convert.label2Dummy(atlas_volume)>0, [],  nRoi); 

%% 2 Wrangle imaging_data

    imaging_shape = size(imaging_volume);
    imaging_volume = reshape(imaging_volume, prod(imaging_shape(1:3)), imaging_shape(4)); % Reshape to 2D   


%% 3 Define sliding window indices

    % define sliding window  
    IX=false(size(imaging_volume,2),1); IX(1:W)=1; % Define and expand window
    IX=cell2mat(arrayfun(@(x) circshift(IX, x), [0:sum(IX):size(IX,1)-(sum(IX))], 'uni',0));
    % Number of windows
    nWin=size(IX,2);
    
%% Calculate FC matrix for each window

    % Pairwise mean timecourse FC for each window
    FC=zeros(nRoi,nRoi,size(IX,2));
    for win = 1:nWin
        FC(:,:,win) = corr(imaging_volume(:,IX(:,win))'*atlas_volume ./ sum(atlas_volume)); % pairwise FC
    end

%% Wrangle FC matrices

    % Retain pos. weightings
    FC = FC .* (FC > 0);

    % Fisher transform r to z-score
    FC=c3nl_fisher(c3nl_symAdj(FC, 'mean'), 'r2z'); 

    % Replace inf's on diagonal with 1
    FC(logical(repmat(eye(size(FC,1)),1,1,nWin)))=1;
    % Replace NaN with 0
    FC(isnan(FC))=0;
    
    
%% Calc dynamic modularity 


    % preallocate
    A = cell(1, nWin);
    B = spalloc(nRoi * nWin, nRoi * nWin,(nRoi + nWin) * nRoi* nWin);
    twomu = 0;
    
    Q = zeros(nRep, 1);


    %--- null model -------------------------------------------------------
    for win = 1 : nWin
        A{win} = FC(:,:,win); % extract FC for a window
        k = sum(A{win});                             % node degree
        twom = sum(k);                               % mean network degree
        twomu = twomu + twom;                        % increment
        indx = [1:nRoi] + (win-1)*nRoi;              % find indices
        B(indx,indx) = A{win} - gamma * [k'*k]/twom; % fill B matrix
    end
    twomu = twomu + 2*omega* nRoi*(nWin-1);
    
    B = B + omega/2*spdiags(ones(nRoi*nWin,2),[-nRoi, nRoi], nRoi*nWin, nRoi*nWin);
    B = B + omega*spdiags(ones(nRoi*nWin,2),[-2*nRoi, 2*nRoi], nRoi*nWin, nRoi*nWin);
    
    for rep = 1 : nRep
        [S,Qtmp] = genlouvain(B);
        Qtmp = Qtmp / twomu;
        S = reshape(S, nRoi, nWin);
        
        Q(rep) = Qtmp;
        modules(:,:,rep) = S;
    end

end