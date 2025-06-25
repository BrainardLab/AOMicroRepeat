clc;
close all;

%---------------------------------------------------------------------------------------------------------
% load 'green_AOM_LUT_processing.mat'

for i = 1:length(all_trials_unpacked)
    lin_intensity = 10.^(all_trials_unpacked(i,1));
    new_intensity(i,1) = green_AOM_lut(round(lin_intensity*1000)+1,2);
%   log_int(i,1)  =  green_AOM_lut(round(new_intensity(i)*255),3)
end

bitmap_intensity = round(255*new_intensity);
all_trials_unpacked(:,4)=bitmap_intensity ;
%green_AOM_lut(:,4) = bitmap_intensity;
orig_int = green_AOM_lut(:,3);

for i = 1:length(all_trials_unpacked)
   log_idx(i) = find(green_AOM_lut(:,4)==bitmap_intensity(i),1);
   log_int_final(i) = orig_int(log_idx(i));
end
log_int_final = log_int_final';

all_trials_unpacked(:,5) = log_int_final;
new_array_int = sortrows(all_trials_unpacked,5);
new_array_int(:,1) = new_array_int(:,5);
new_array_int =  new_array_int(:,1:4);
save LUTcorrected_vector new_array_int
