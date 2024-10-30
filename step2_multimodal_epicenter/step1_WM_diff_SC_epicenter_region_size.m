%% find the epicenter of WM with HCP WM co-variance matrix in BD and SCH
path = pwd;

cd ..
load('region_size_JHU.mat')
load('node_name.mat')
load('Results_signal.mat')
diff_results = results;

path_data = [pwd '\public_data'];

cd('step1_multimodal_corr')
rand_path = [pwd '\Burt\data_rand'];

Type = {'volume'};
Group = {'BD','SCH'};
Group_name = {'BD','MDD','SCH'};

roi_num = 48;
rand_num = 10000;
p_thr = 0.05;

results = struct;
results(1) = [];
results_abs = results;


%% cal
for itype = 1:length(Type)
    cd(path_data)
    load('unrelated_volume_JHU.mat','SC')

    for igroup = 1:length(Group)
        group_index = 0;
        for iname = 1:length(Group_name)
            if isequal(Group_name{iname},Group{igroup}) ...
                    && group_index == 0
                group_index = iname;
            end
        end

        results(end+1).type = Type{itype};
        results(end).group = Group{igroup};
        results(end).Z = diff_results.([Type{itype} '_Z'])(:,group_index);
        results(end).Net = SC;

        results_abs(end+1).type = Type{itype};
        results_abs(end).group = Group{igroup};
        results_abs(end).Z = diff_results.([Type{itype} '_Z'])(:,group_index);
        results_abs(end).Net = abs(SC);

        clear surrogate_maps
        cd([rand_path '\' Type{itype} '_' Group{igroup}])
        load([Type{itype} '_' Group{igroup} '_burt_rand.mat'])
        surrogate_maps = surrogate_maps(1:rand_num,:);

        for iroi = 1:roi_num
            disp(['Now calculating the roi(' num2str(iroi) '/' num2str(roi_num) ...
                ') in ' Type{itype} '_' Group{igroup} '  |' datestr(clock)])

            roi_index = setdiff(1:roi_num,iroi);

            results(end).rho(iroi,1) = ...
                partialcorr(results(end).Z(roi_index),SC(roi_index,iroi),region_size(roi_index),'type','Spearman');
            results_abs(end).rho_abs(iroi,1) = ...
                partialcorr(results_abs(end).Z(roi_index),abs(SC(roi_index,iroi)),region_size(roi_index),'type','Spearman');

            results(end).rho_rand(iroi,:) = ...
                partialcorr(surrogate_maps(:,roi_index)',SC(roi_index,iroi),region_size(roi_index),'type','Spearman');
            results_abs(end).rho_abs_rand(iroi,:) = ...
                partialcorr(surrogate_maps(:,roi_index)',abs(SC(roi_index,iroi)),region_size(roi_index),'type','Spearman');
        end

        results(end).rho_p = (sum(repmat(results(end).rho,1,rand_num) >= results(end).rho_rand,2)+1)/(rand_num+1);
        results_abs(end).rho_abs_p = (sum(repmat(results_abs(end).rho_abs,1,rand_num) >= results_abs(end).rho_abs_rand,2)+1)/(rand_num+1);

        [~,I] = sort(results(end).rho,'ascend');
        results(end).rho_sort = [node_name(I),num2cell(results(end).rho(I)),num2cell(results(end).rho_p(I))];
        [~,I] = sort(results_abs(end).rho_abs,'ascend');
        results_abs(end).rho_abs_sort = [node_name(I),num2cell(results_abs(end).rho_abs(I)),num2cell(results_abs(end).rho_abs_p(I))];
    end
end


%% save
cd(path)
save('Results_epicenter_diff_SC_region_size.mat','results','results_abs')