%% Function to calculate mean and differences
%
% This is useful for creating Bland-Altman plots
function [mean_val, diff_val] = calculate_bland_altman(x1, x2)
    mean_val = (x1 + x2) / 2;
    diff_val = x1 - x2;
end