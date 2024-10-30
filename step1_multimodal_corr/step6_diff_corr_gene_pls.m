%% corr to gene
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


%% gene expression
cd(path_data)
data_gene = readmatrix('gene_JHU.csv');
data_gene(:,1) = [];

nanroi = isnan(data_gene(:,1));
data_gene = data_gene(~nanroi,:);


%% corr
for itype = 1:length(Type)
    % signal diff
    for igroup = 1:length(group)
        results(end+1).type = Type{itype};
        results(end).group = group{igroup};
        results(end).Z = results_signal.([Type{itype} '_Z'])(~nanroi,igroup);
        results(end).signal = data_gene;

        [~,~,XS,~,~,PCTVAR,~,stats] = plsregress(results(end).signal,results(end).Z,1);
        results(end).XS = XS;
        results(end).PCTVAR = PCTVAR;
        results(end).W = stats.W;

        results(end).PCTVAR_rand = zeros(rand_time,1);
        cd(rand_path)
        cd([Type{itype} '_' group{igroup}(1:end-1)])
        clear surrogate_maps
        load([Type{itype} '_' group{igroup}(1:end-1) '_burt_rand.mat'])
        for irand = 1:rand_time
            if mod(irand,floor(rand_time/10)) == 0
                disp(['Now processing the data in ' Type{itype} '_signal_' group{igroup} ...
                    num2str(irand) '\' num2str(rand_time) '  |' datestr(clock)])
            end

            [~,~,~,~,~,PCTVAR,~,stats] = plsregress(results(end).signal,surrogate_maps(irand,~nanroi)',1);
            results(end).PCTVAR_rand(irand) = PCTVAR(2);
        end
        results(end).p = (sum(abs(results(end).PCTVAR(2)) <= abs(results(end).PCTVAR_rand)) + 1)/(rand_time + 1);
    end
end


%% save
cd(path)
save('Results_diff_corr_gene_all_pls.mat','results')