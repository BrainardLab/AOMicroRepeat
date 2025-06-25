% Clean the workspace
clc;
close all;
clear all;
%-----------------------------------------------------------------------------
% Get the file path
path = 'W:\Data\11125\20231103\AO_Psychophysics\QUEST\43\group2';
% (specify the folder path in this format to avoid mixing of trials with
% different stimulus sizes) Also check if all the subfolders are added to
% the current path.

trial_videos = dir(path);
num_trial_videos = length(trial_videos(3:end));

%-----------------------------------------------------------------------------
%Extract the matfiles from all the trials of a single stimulus size (8 or
%43)
list_mat_file_path={};
mat_file={}; 
all_trials = {};
trial_intervals = [44,88,132,176];

for i = 3: size(trial_videos) 
       i
       list_mat_file_path{i} = fullfile(path,trial_videos(i).name);
       cd(list_mat_file_path{1,i})
       curr_path = cd(list_mat_file_path{1,i});
       mat_file{i,1} = (dir(fullfile(curr_path,'*.mat')));
       Trial = (mat_file{i});
       Trial_unpack = Trial(2,1).name %sometines (1,1)
       load(Trial_unpack)
       all_trials{i,1}=trial_matrix;
       all_trials{i,2}=response_matrix;
       all_trials{i,3}=theThreshold;
          
end

%----------------------------------------------------------------------------
% Unpack Matfiles here

all_trials_unpacked = nan(88,3); % total number of trials for QUEST for the current experiemnt(44x4)

% Unpacking Intensities here
all_trials_unpacked(1:44,1)= all_trials{3,1};
all_trials_unpacked(45:88,1)= all_trials{4,1};
% all_trials_unpacked(89:132,1)= all_trials{5,1};
% all_trials_unpacked(133:176,1)= all_trials{6,1};

% Unpacking Responses here
all_trials_unpacked(1:44,2)= all_trials{3,2};
all_trials_unpacked(45:88,2)= all_trials{4,2};
% all_trials_unpacked(89:132,2)= all_trials{5,2};
% all_trials_unpacked(133:176,2)= all_trials{6,2};

% Unpacking Thresholds here
all_trials_unpacked(1:44,3)= all_trials{3,3};
all_trials_unpacked(45:88,3)= all_trials{4,3};
% all_trials_unpacked(89:132,3)= all_trials{5,3};
% all_trials_unpacked(133:176,3)= all_trials{6,3};

%-----------------------------------------------------------------------------
% Save the all trials mat file in the original path
cd(path)
save all_trials_unpacked.mat
save 11108_QUEST_8 all_trials_unpacked
%-----------------------------------------------------------------------------