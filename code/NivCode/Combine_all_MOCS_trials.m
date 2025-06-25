% Clean the workspace
clc;
close all;
clear all;
%-----------------------------------------------------------------------------
% Get the file path
path = 'W:\Data\11125\20231103\AO_Psychophysics\MOCS\8\group2';
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


for i = 3: size(trial_videos) 
       i
       list_mat_file_path{i} = fullfile(path,trial_videos(i).name);
       cd(list_mat_file_path{1,i})
       curr_path = cd(list_mat_file_path{1,i});
       mat_file{i,1} = (dir(fullfile(curr_path,'*.mat')));
       Trial = (mat_file{i});
       Trial_unpack = Trial(2,1).name;
       load(Trial_unpack)
       all_trials{i,1}=trial_vector;
       all_trials{i,2}=response_vector;   
end

%----------------------------------------------------------------------------
% Unpack Matfiles here

all_trials_unpacked = nan(length(all_trials{3,1}),2); % total number of trials for MOCS in the current experiemnt(100 or 90)

% Unpacking Intensities here
if length(all_trials{3,1}) == 90

  all_trials_unpacked(1:90,1)= all_trials{3,1};
  all_trials_unpacked(91:180,1)= all_trials{4,1};
%   all_trials_unpacked(181:270,1)= all_trials{5,1};
%   all_trials_unpacked(271:360,1)= all_trials{6,1};

% Unpacking Responses here
  all_trials_unpacked(1:90,2)= all_trials{3,2};
  all_trials_unpacked(91:180,2)= all_trials{4,2};
%   all_trials_unpacked(181:270,2)= all_trials{5,2};
%   all_trials_unpacked(271:360,2)= all_trials{6,2};
  
else 
   all_trials_unpacked(1:100,1)= all_trials{3,1};
  all_trials_unpacked(101:200,1)= all_trials{4,1};
%   all_trials_unpacked(201:300,1)= all_trials{5,1};
%   all_trials_unpacked(301:400,1)= all_trials{6,1};

% Unpacking Responses here
  all_trials_unpacked(1:100,2)= all_trials{3,2};
  all_trials_unpacked(101:200,2)= all_trials{4,2};
%   all_trials_unpacked(201:300,2)= all_trials{5,2};
%   all_trials_unpacked(301:400,2)= all_trials{6,2};
 
end

   
%-----------------------------------------------------------------------------
% Save the all trials mat file in the original path
cd(path)
save all_trials_unpacked.mat %(note:intensities are not sorted here by acsending/decending order)
%-----------------------------------------------------------------------------