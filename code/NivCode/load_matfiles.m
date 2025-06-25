% Clean the workspace
clc;
close all;
clear all;
%-----------------------------------------------------------------------------
% Get the file path
path = 'W:\Data\11119\20231012\AO_Psychophysics\MOCS\8';
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

       list_mat_file_path{i} = fullfile(path,trial_videos(i).name);
       curr_path = cd(list_mat_file_path{1,i});
       mat_file{i,1} = (dir(fullfile(curr_path,'*.mat')));
       Trial1 = (mat_file{i});
       Trial_1 = Trial1(2,1).name;
       load(Trial_1)
       all_trials{i,1}=trial_vector;
       all_trials{i,2}=response_vector;
%        all_trials{i,3}=theThreshold;
          
end

%----------------------------------------------------------------------------
% Unpack Matfiles here

all_trials_unpacked = nan(176,3); % total number of trials for QUEST for the current experiemnt(44x4)

% Unpacking Intensities here
all_trials_unpacked(1:44,1)= all_trials{3,1};
all_trials_unpacked(45:88,1)= all_trials{4,1};
all_trials_unpacked(89:132,1)= all_trials{5,1};
all_trials_unpacked(133:176,1)= all_trials{6,1};

% Unpacking Responses here
all_trials_unpacked(1:44,2)= all_trials{3,2};
all_trials_unpacked(45:88,2)= all_trials{4,2};
all_trials_unpacked(89:132,2)= all_trials{5,2};
all_trials_unpacked(133:176,2)= all_trials{6,2};

% Unpacking Thresholds here
all_trials_unpacked(1:44,3)= all_trials{3,3};
all_trials_unpacked(45:88,3)= all_trials{4,3};
all_trials_unpacked(89:132,3)= all_trials{5,3};
all_trials_unpacked(133:176,3)= all_trials{6,3};

%-----------------------------------------------------------------------------
% Save the all trials mat file in the original path
cd(path)
save all_trials_unpacked.mat
%-----------------------------------------------------------------------------