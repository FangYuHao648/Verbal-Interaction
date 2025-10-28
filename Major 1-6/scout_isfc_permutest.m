function critical_t_value = scout_isfc_permutest(data_1,data_2)
% input:
%           data_1: task condition data
%           data_2: resting state data
%  data_1 and data_2: 62 rois x 62 rois x 29 subject pairs 
% output:
%          critical_t_value: 临界t值

% permutation
n_permutations = 500;
max_cluster_t_stats = zeros(n_permutations, 1);
mergedata=cat(3,data_1,data_2);
nv = size(data_1,1); % number of vertices/ channels / rois
nsbj = size(data_1,3);

idx = 1:1:2*nsbj;
for perm = 1:n_permutations
    % 随机打乱标签
    perm_indices = idx(randperm(numel(idx),nsbj));
    perm_data_1 = mergedata(:, :, perm_indices);
    perm_data_2 = mergedata(:, :, idx(~ismember(idx,perm_indices)));

    % 重新计算 t 统计量
    perm_t_values = zeros(nv);
    for i = 1:nv 
        for j = 1:nv 
            % 提取连接强度
            perm_edge_1 = squeeze(perm_data_1(i, j, :));
            perm_edge_2 = squeeze(perm_data_2(i, j, :));

            % 计算 t 统计量
            [~, ~, ~, stats] = ttest(perm_edge_1, perm_edge_2);
            perm_t_values(i, j) = stats.tstat;
        end
    end

    % 筛选显著边缘
    perm_significant_edges_positive = perm_t_values > 1.96;
    perm_significant_edges_negative = perm_t_values < -1.96;

    perm_bipartite_adj_positive = zeros(2*nv);
    perm_bipartite_adj_positive(1:nv, nv+1:2*nv) = perm_significant_edges_positive;
    perm_bipartite_adj_positive(nv+1:2*nv, 1:nv) = perm_significant_edges_positive';

    perm_bipartite_adj_negative = zeros(2*nv);
    perm_bipartite_adj_negative(1:nv, nv+1:2*nv) = perm_significant_edges_negative;
    perm_bipartite_adj_negative(nv+1:2*nv, 1:nv) = perm_significant_edges_negative';

    G_perm_positive = graph(perm_bipartite_adj_positive);
    G_perm_negative = graph(perm_bipartite_adj_negative);

    % 查找所有连通分量
    bins_perm_positive = conncomp(G_perm_positive);
    num_bins_perm_positive = max(bins_perm_positive);
    bins_perm_negative = conncomp(G_perm_negative);
    num_bins_perm_negative = max(bins_perm_negative);

    % 计算每个连通分量的 t 统计量总和
    max_t_stat_sum = 0;
    for k = 1:num_bins_perm_positive
        cluster_nodes_positive = find(bins_perm_positive == k);
        t_stat_sum = sum(sum(abs(perm_t_values(cluster_nodes_positive(cluster_nodes_positive <= nv), cluster_nodes_positive(cluster_nodes_positive > nv) - nv))));
        if t_stat_sum > max_t_stat_sum
            max_t_stat_sum = t_stat_sum;
        end
    end
    for k = 1:num_bins_perm_negative
        cluster_nodes_negative = find(bins_perm_negative == k);
        t_stat_sum = sum(sum(abs(perm_t_values(cluster_nodes_negative(cluster_nodes_negative <= nv), cluster_nodes_negative(cluster_nodes_negative > nv) - nv))));
        if t_stat_sum > max_t_stat_sum
            max_t_stat_sum = t_stat_sum;
        end
    end
    max_cluster_t_stats(perm) = max_t_stat_sum;
end

% 计算临界 t 值
critical_t_value = prctile(max_cluster_t_stats, 95);

end