%%% Confine/bin MOCS 43 results to -3.

intensities =  all_trials_unpacked(:,1);

for i  = 1: length(intensities)
    if intensities(i) <-3 && intensities(i) >=-3.5
        intensities(i) = -3.5;
    else
        %do nothing
    end
end

all_trials_unpacked = [intensities,all_trials_unpacked(:,2),all_trials_unpacked(:,3)];
save 11002_S2_SS8 all_trials_unpacked

