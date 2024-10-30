%% corr to WM volume
%% para
path = pwd;
rand_path = [path '\Burt\data_rand'];
rand_time = 10000;

cd ..
load('Results_signal.mat')
results_signal = results;

cd('public_data')
path_data = pwd;

Type = {'volume'};
group = {'BD ','MDD ','SCH '};
p_thr = 0.05;

results = struct;
results(1) = [];

%% corr
for itype = 1:length(Type)
    % signal data
    cd(path_data)
    data_signal = load(['unrelated_' Type{itype} '_JHU.mat']);
    data_signal = data_signal.([Type{itype} '_signal']);
    data_signal = mean(data_signal,1)';

    % signal diff
    results(end+1).type = Type{itype};
    results(end).group = cat(2,group{:});
    results(end).Z = results_signal.([Type{itype} '_Z']);
    results(end).signal = data_signal;
    results(end).rho = corr(results(end).Z,results(end).signal,'type','Spearman');

    results(end).rho_rand = zeros(size(results(end).rho,1),rand_time);
    for igroup = 1:length(group)
        cd(rand_path)
        cd([Type{itype} '_' group{igroup}(1:end-1)])

        clear surrogate_maps
        load([Type{itype} '_' group{igroup}(1:end-1) '_burt_rand.mat'])
        results(end).rho_rand(igroup,:) = corr(surrogate_maps',results(end).signal,'type','Spearman');
    end

    results(end).p = (sum(repmat(abs(results(end).rho),1,rand_time) <= abs(results(end).rho_rand),2) + 1)/(rand_time + 1);
    results(end).sig = cat(2,group{results(end).p <= p_thr});
end


%% save
cd(path)
save('Results_diff_corr_volume.mat','results')