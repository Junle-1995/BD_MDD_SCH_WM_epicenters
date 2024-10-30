%% find the epicenter of WM using Lp with HCP T1/T2 co-variance matrix in BD and SCH
path = pwd;

load('Results_epicenter_diff_SC_T1_T2_region_size.mat')
results_con = results;
results_con_abs = results_abs;

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


%% Net
cd(path_data)
load('unrelated_T1_T2_ratio_JHU.mat','SC')

SC_ori = SC;
[SC, ~] = gretna_batch_lp_bin_all(SC_ori, 0.09, 0.3, 0.02, 'pos', 's');
[SC_abs, ~] = gretna_batch_lp_bin_all(abs(SC_ori), 0.09, 0.3, 0.02, 'pos', 's');


%% cal
for itype = 1:length(Type)
    for igroup = 1:length(Group)
        % find the index of connectivity corr results
        i_index = [];
        for i = 1:length(results_con)
            if isequal(results_con(i).type,Type{itype}) && ...
                    isequal(results_con(i).group,Group{igroup})
                if ~isempty(i_index)
                    error('?')
                else
                    i_index = i;
                end
            end
        end

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
        results_abs(end).Net = SC_abs;

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
                partialcorr(results_abs(end).Z(roi_index),SC_abs(roi_index,iroi),region_size(roi_index),'type','Spearman');

            results(end).rho_rand(iroi,:) = ...
                partialcorr(surrogate_maps(:,roi_index)',SC(roi_index,iroi),region_size(roi_index),'type','Spearman');
            results_abs(end).rho_abs_rand(iroi,:) = ...
                partialcorr(surrogate_maps(:,roi_index)',SC_abs(roi_index,iroi),region_size(roi_index),'type','Spearman');
        end

        results(end).rho_p = (sum(repmat(results(end).rho,1,rand_num) <= results(end).rho_rand,2)+1)/(rand_num+1);
        results_abs(end).rho_abs_p = (sum(repmat(results_abs(end).rho_abs,1,rand_num) <= results_abs(end).rho_abs_rand,2)+1)/(rand_num+1);

        [~,I] = sort(results_con(i_index).rho,'ascend');
        results(end).rho_sort = [node_name(I),num2cell(results(end).rho(I)),num2cell(results(end).rho_p(I))];
        [~,I] = sort(results_con_abs(i_index).rho_abs,'ascend');
        results_abs(end).rho_abs_sort = [node_name(I),num2cell(results_abs(end).rho_abs(I)),num2cell(results_abs(end).rho_abs_p(I))];
    end
end


%% save
cd(path)
save('Results_epicenter_diff_SC_T1_T2_Lp_region_size.mat','results','results_abs')