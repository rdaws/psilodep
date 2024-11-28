function [Q_norm, FC, TC] = calculate_static_modularity(imaging_volume, atlas_volume, nP)
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
% Psilodep :- Calculate modularity (Q)
% 
% Requires 
%    third_party folder and subfolders to be in path 
%
%
%   Args:
%       imagine_volume - 4D matrix of fMRI imaging data
%       atlas_volume   - 3D matrix of unique values denoting Regions of
%                        Interest(ROI)
%       nP             - Number of permutations to run when estimating
%                        modularity
%   Returns:
%       Q_norm         - Single value representing normalised Q
%       FC             - 2D functional connectivity matrix
%       TC             - 3D matrix of ROI timecourses
%
%
% Steps
% 1 - Wrangle atlas_volume to be 2D
% 2 - Wrangle imaging_volume to be 2D
% 3 - Extract mean timecourses for each atlas region
% 4 - Define FC matrix
% 5 - Estimate Modularity (Q)
% 6 - Normalise Q
%
%
% Adapted from Karolina Finc et al 2020, Nat Comms.
%

% Number of permutations
if (~exist('nP', 'var'))
    nP = 100;
end
    


%% 1 Wrangle atlas

    % Find number of ROIs
    nRoi=numel(unique(atlas_volume(atlas_volume>0)));
    % Reshape atlas to 2D
    atlas_volume=reshape(label2Dummy(atlas_volume)>0, [],  nRoi); 

%% 2 Wrangle imaging_data

    imaging_shape = size(imaging_volume);
    imaging_volume = reshape(imaging_volume, prod(imaging_shape(1:3)), imaging_shape(4)); % Reshape to 2D     


%% 3 Extract mean timecourses for each atlas region

    TC = imaging_volume' * atlas_volume ./ sum(atlas_volume); % Extract mean timecourse, per roi


%% 4 define FC matrix with Pearson correlation and filtering

    % Define functional connectivity with Pearson correlation
    FC = corr(TC);

    % Remove negative correlations
    FC = FC .* (FC > 0);

    % Fisher transform r to z-scores
    FC=c3nl_fisher(c3nl_symAdj(FC, 'mean'), 'r2z');

    % Replace inf's on diagonal with 1
    FC(logical(eye(size(FC))))=1;

    % Replace NaN with 0
    FC(isnan(FC))=0;


%% 5 Estimate Modularity (Q) for the real and shuffled FC matrices nP times

    Q = zeros(nP, 1);
    Q_null = zeros(nP, 1);

    % Run Q estimation on real and shuffled networks nP times
    for ii = 1:nP
        [~,Q(ii)]=community_louvain(FC,1);
        [~,Q_null(ii)]=community_louvain(randmio_und(FC,1),1);
    end


%% 6 Normalise Q

    % Divide the max real Q value by the mean of the shuffled Q values
    Q_norm = max(Q) / mean(Q_null);



end
