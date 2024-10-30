function stepX_WM_epi_SC_T1_T2_corr_GM_diff_Lp_region_size(GM_atlas,GM_roi_num,s_min,s_max)
%% corr the GM volume diff Z map with Lp of SC T1_T2 to the WM epicenter
path = pwd;

cd ..
load(['Results_signal_GM_' GM_atlas '.mat'])
diff_results = results;

cd('step2_multimodal_epicenter')
load('Results_epicenter_diff_SC_T1_T2_Lp_region_size.mat')
lp_results = results;

cd ..
path_data = [pwd '\public_data'];
load('node_name.mat')

% GM_atlas = 'Sc400';
load(['region_size_' GM_atlas '.mat'])
% GM_roi_num = 400;
JHU_roi_num = 48;
Group_name = {'BD','MDD','SCH'};

cd(path)
rand_path = [pwd '\Burt\data_rand'];

rand_num = 10000;
p_thr = 0.05;

results = struct;


%% cal
for i = 1:length(lp_results)
    cd(path_data)
    load(['unrelated_SC_T1_T2_' GM_atlas '_JHU.mat'],'SC')

    cd ..
    [Net, ~] = gretna_batch_lp_bin_all(SC, s_min, s_max, 0.02, 'pos', 's');

    results(i).type = lp_results(i).type;
    results(i).group = lp_results(i).group;
    results(i).epi_index = find(lp_results(i).rho_p <= p_thr/JHU_roi_num);
    results(i).epi_name = node_name(results(i).epi_index);

    if isempty(results(i).epi_index)
        continue,
    end

    group_index = 0;
    for igroup = 1:length(Group_name)
        if isequal(Group_name{igroup},results(i).group) ...
                && group_index == 0
            group_index = igroup;
        end
    end

    data_Z = diff_results.([results(i).type '_Z'])(:,group_index);
    results(i).data_T = data_Z;
    results(i).region_size = region_size;

    clear surrogate_maps
    cd([rand_path '\' GM_atlas '_' results(i).type '_' results(i).group])
    load([GM_atlas '_' results(i).type '_' results(i).group '_burt_rand.mat'])
    surrogate_maps = surrogate_maps(1:rand_num,:);

    for iroi = 1:length(results(i).epi_index)
        epi = results(i).epi_index(iroi);
        results(i).data_net(:,iroi) = Net(1:GM_roi_num,epi+GM_roi_num);

        results(i).rho(iroi,1) = ...
            partialcorr(data_Z,Net(1:GM_roi_num,epi+GM_roi_num),region_size,'type','Spearman');

        results(i).rho_rand(iroi,:) = ...
            partialcorr(surrogate_maps',Net(1:GM_roi_num,epi+GM_roi_num),region_size,'type','Spearman');
    end

    results(i).rho_p = (sum(repmat(results(i).rho,1,rand_num) <= results(i).rho_rand,2)+1)/(rand_num+1);
    results(i).rho_p_Bon = p_thr/length(results(i).rho_p);
    
    results(i).rho_pattern = [results(i).epi_name,num2cell(results(i).rho),num2cell(results(i).rho_p)];
end


%% save
cd(path)
save(['Results_WM_epicenter_SC_T1_T2_' GM_atlas '_diff_Lp_region_size.mat'],'results')

return