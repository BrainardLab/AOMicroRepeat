
% Load a QUEST data file and run this to find maximum trial intensity with
% incorrect response. Intensity is raw log10, not LUT corrected.
index = find(response_matrix == 0 & ~isnan(theThreshold) );
max(trial_matrix(index))

% Load a MOCS data file and run this to find maximum trial intensity with
% incorrect response.  Intensity is raw log10, not LUT corrected.
index = find(response_vector);
max(trial_vector(index))