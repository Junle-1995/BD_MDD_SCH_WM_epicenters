%% corr to fun TC std
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


%% cal the BOLD TC std
cd(path_data)
load('unrelated_fun_TC_std.mat')
data_tc = mean(data_tc,2);


%% corr
for itype = 1:length(Type)
    % signal diff
    results(end+1).type = Type{itype};
    results(end).group = cat(2,group{:});
    results(end).Z = results_signal.([Type{itype} '_Z']);
    results(end).signal = data_tc;
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
    results(end).p_fdr = gretna_FDR(results(end).p,p_thr);
    results(end).sig = cat(2,group{results(end).p <= p_thr});
end


%% save
cd(path)
save('Results_diff_corr_fun_tc_std.mat','results')