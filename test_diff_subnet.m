%% test the Z-value of each subnet
path = pwd;

load('Results_signal.mat')
load('subnet.mat')

Type = {'volume'};
Group = {'BD','MDD','SCH'};

rand_time = 10000;
p_thr = 0.05;

results_subnet = struct;
results_subnet(1) = [];


%% cal
for igroup = 1:length(Group)
    for itype = 1:length(Type)
        Z = results.([Type{itype} '_Z'])(:,igroup);

        results_subnet(end+1).Group = Group{igroup};
        results_subnet(end).Type = Type{itype};

        for inet = 1:length(subnet_name)
            results_subnet(end).subnet_Z_cell{inet,1} = abs(Z(subnet_index{inet}));
            results_subnet(end).subnet_Z(inet,1) = mean(abs(Z(subnet_index{inet})));
        end

        for irand = 1:rand_time
            rand_order = randperm(length(Z));
            while isequal(rand_order,1:length(Z))
                rand_order = randperm(length(Z));
            end
            Z_rand = Z(rand_order);

            for inet = 1:length(subnet_name)
                results_subnet(end).subnet_Z_rand(inet,irand) = mean(abs(Z_rand(subnet_index{inet})));
            end
        end

        results_subnet(end).subnet_P = ...
            (sum(repmat(results_subnet(end).subnet_Z,1,rand_time) <= results_subnet(end).subnet_Z_rand,2)+1)/(rand_time+1);

        results_subnet(end).subnet_sig = subnet_name(results_subnet(end).subnet_P <= p_thr/length(results_subnet(end).subnet_P));
    end
end

save('Results_diff_subnet_test.mat','results_subnet')