%% test the subnet to epicenter
path = pwd;

cd ..
load('subnet.mat')

Type = {'FC','SC','SC_T1_T2','TC','FC_Lp','SC_Lp','SC_T1_T2_Lp','TC_Lp'};

p_thr = 0.05;
rand_time = 10000;

subnet_results = struct;
subnet_results(1) = [];

cd(path)


%% cal
for itype = 1:length(Type)
    load(['Results_epicenter_diff_' Type{itype} '_region_size.mat'])

    for igroup = 1:length(results)
        epi_index = find(results(igroup).rho_p <= p_thr/length(results(igroup).rho_p));
        if ~isempty(epi_index)
            subnet_results(end+1).Type = Type{itype};
            subnet_results(end).feature = results(igroup).type;
            subnet_results(end).group = results(igroup).group;
            subnet_results(end).epi_index = epi_index;
            subnet_results(end).subnet_name = subnet_name;

            for inet = 1:length(subnet_index)
                subnet_results(end).subnet_num(inet,1) = ...
                    length(intersect(subnet_index{inet},epi_index));

                subnet_results(end).subnet_percent(inet,1) = ...
                    subnet_results(end).subnet_num(inet,1)/length(subnet_index{inet});
            end

            for irand = 1:rand_time
                rand_order = randperm(length(results(igroup).rho_p));
                while isequal(rand_order,1:length(results(igroup).rho_p))
                    rand_order = randperm(length(results(igroup).rho_p));
                end
                epi_index_rand = rand_order(epi_index);

                for inet = 1:length(subnet_index)
                    subnet_results(end).subnet_num_rand(inet,irand) = ...
                        length(intersect(subnet_index{inet},epi_index_rand));
                end
            end

            subnet_results(end).subnet_num_p = ...
                (sum(repmat(subnet_results(end).subnet_num,1,rand_time) <= subnet_results(end).subnet_num_rand,2)+1)/(rand_time+1);

            subnet_results(end).subnet_sig = subnet_name(subnet_results(end).subnet_num_p <= p_thr/length(subnet_results(end).subnet_num_p))';
        end
    end
end

save('Results_epicenter_diff_subnet_region_size.mat','subnet_results')