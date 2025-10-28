function [significant_edges_positive, significant_edges_negative,significant_cluster_matrices,max_pos,pos_idx,max_neg,neg_idx] = scout_significant_edges_test(t_values, data_1, data_2)
% input:
%           data_1: task condition data
%           data_2: resting state data
%           t_value: t_value of paired edges between data_1 and data_2
%  data_1 and data_2: 62 rois x 62 rois x 29 subject pairs


nv = size(t_values,1);

% 筛选显著边缘
significant_edges_positive = t_values > 1.96;
significant_edges_negative = t_values < -1.96; % 这里有t对负值的判断，line46 adding correct value是否正确？

% 创建二分图的邻接矩阵，下方的正负是根据t值的正负来区分
bipartite_adj_positive = zeros(2*nv);
bipartite_adj_positive(1:nv, nv+1:2*nv) = significant_edges_positive;  % right-top block
bipartite_adj_positive(nv+1:2*nv, 1:nv) = significant_edges_positive'; % left-bottom block

bipartite_adj_negative = zeros(2*nv);
bipartite_adj_negative(1:nv, nv+1:2*nv) = significant_edges_negative;  % right-top block
bipartite_adj_negative(nv+1:2*nv, 1:nv) = significant_edges_negative'; % left-bottom block
% 创建无向二分图
G_positive = graph(bipartite_adj_positive);
G_negative = graph(bipartite_adj_negative);

% permutation
critical_t_value = scout_isfc_permutest(data_1,data_2);


% 初始化显著簇的矩阵
significant_cluster_matrices = {};

% 计算实际数据的簇 t 统计量
actual_clusters_positive = conncomp(G_positive); % 返回经过第几个电极的强连通分图SCC的个数
num_actual_clusters_positive = max(actual_clusters_positive);
significant_clusters_positive = [];

actual_clusters_negative = conncomp(G_negative);
num_actual_clusters_negative = max(actual_clusters_negative);
significant_clusters_negative = [];

for k = 1:num_actual_clusters_positive
    % 返回SCC=k时的电极标号，但范围是1：2*nv
    cluster_nodes_positive = find(actual_clusters_positive == k);

    % t_values nvxnv
    t_stat_sum = sum(sum(abs(t_values(cluster_nodes_positive(cluster_nodes_positive <= nv), cluster_nodes_positive(cluster_nodes_positive > nv) - nv))));
    positive_t_sum(1,k)=t_stat_sum;
    if t_stat_sum > critical_t_value
        significant_clusters_positive = [significant_clusters_positive, {cluster_nodes_positive}];
        % 构建显著簇的矩阵
        cluster_matrix_positive = zeros(2*nv, 2*nv);
        for node1 = cluster_nodes_positive
            if node1 <= nv
                for node2 = cluster_nodes_positive
                    if node2 > nv
                        if significant_edges_positive(node1, node2 - nv)
                            cluster_matrix_positive(node1, node2) = t_values(node1, node2 - nv);
                            cluster_matrix_positive(node2, node1) = t_values(node1, node2 - nv);
                        end
                    end
                end
            end
        end
        significant_cluster_matrices{end + 1} = cluster_matrix_positive;
    end
end
% max_pos: maximum positive t value sum, pos_idx: SCC的实际个数
[max_pos,pos_idx]=max(positive_t_sum);

for n = 1:num_actual_clusters_negative
    cluster_nodes_negative = find(actual_clusters_negative == n);
    t_stat_sum = sum(sum(abs(t_values(cluster_nodes_negative(cluster_nodes_negative <= nv), cluster_nodes_negative(cluster_nodes_negative > nv) - nv))));
    negative_t_sum(1,n)=t_stat_sum;
    if t_stat_sum > critical_t_value
        significant_clusters_negative = [significant_clusters_negative, {cluster_nodes_negative}];
        % 构建显著簇的矩阵
        cluster_matrix_negative = zeros(2*nv);
        for node1 = cluster_nodes_negative
            if node1 <= nv
                for node2 = cluster_nodes_negative
                    if node2 > nv
                        if significant_edges_negative(node1, node2 - nv)
                            cluster_matrix_negative(node1, node2) = t_values(node1, node2 - nv);
                            cluster_matrix_negative(node2, node1) = t_values(node1, node2 - nv);
                        end
                    end
                end
            end
        end
        significant_cluster_matrices{end + 1} = cluster_matrix_negative;
    end
end
[max_neg,neg_idx]=max(negative_t_sum);


end