%% corr to neurotransmitter
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


%% transmitter
cd(path_data)
load('neuro_trans_JHU.mat');
neuro_name_all = fieldnames(data_neuro_JHU);

neuro_name = {};
data_neuro = [];
for ineuro = 1:length(neuro_name_all)
    if isempty(strfind(neuro_name_all{ineuro},'CBF'))
        data_neuro(:,end+1) = data_neuro_JHU.(neuro_name_all{ineuro});
        neuro_name{end+1,1} = neuro_name_all{ineuro}(3:end);
    end
end


%% corr
for itype = 1:length(Type)
    % signal diff
    results(end+1).type = Type{itype};
    results(end).group = cat(2,group{:});
    results(end).Z = results_signal.([Type{itype} '_Z']);
    results(end).signal = data_neuro;
    results(end).rho = corr(results(end).Z,results(end).signal,'type','Spearman','rows','pairwise');

    results(end).rho_rand = zeros(size(results(end).rho,1),size(results(end).rho,2),rand_time);
    for igroup = 1:length(group)
        cd(rand_path)
        cd([Type{itype} '_' group{igroup}(1:end-1)])

        clear surrogate_maps
        load([Type{itype} '_' group{igroup}(1:end-1) '_burt_rand.mat'])
        results(end).rho_rand(igroup,:,:) = corr(surrogate_maps',results(end).signal,'type','Spearman','rows','pairwise')';
    end

    results(end).p = (sum(repmat(abs(results(end).rho),1,1,rand_time) <= abs(results(end).rho_rand),3) + 1)/(rand_time + 1);
    
    cd(path)
    cd ..
    p_fdr = gretna_FDR(results(end).p(:),p_thr);
    if ~isempty(p_fdr)
        for igroup = 1:length(group)
            results(end).([group{igroup}(1:end-1) '_sig']) = [neuro_name(results(end).p(igroup,:) <= p_fdr),...
                num2cell(results(end).rho(igroup,results(end).p(igroup,:) <= p_fdr)'),...
                num2cell(results(end).p(igroup,results(end).p(igroup,:) <= p_fdr)')];
        end
    else
        results(end).([group{igroup}(1:end-1) '_sig']) = [];
    end
end


%% save
cd(path)
save('Results_diff_corr_neurotrans.mat','results')