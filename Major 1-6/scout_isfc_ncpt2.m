% reproduce the ncpt as NCPtest.m
% load scout_isfc_all.mat;



F2F_ISFC_ave = mean(F2F_ISFC_all,3);
TWE_ISFC_ave = mean(TWE_ISFC_all,3);
ViC_ISFC_ave = mean(ViC_ISFC_all,3);
bsl_ISFC_ave = mean(bsl_ISFC_all,3);
    figure 
subplot(221)
    imagesc(F2F_ISFC_ave); colorbar; axis image; title(' F2F\_ISFC\_ave (62×62)');
subplot(222) 
    imagesc(TWE_ISFC_ave); colorbar; axis image; title(' TWE\_ISFC\_ave (62×62)');
subplot(223)
    imagesc(ViC_ISFC_ave); colorbar; axis image; title(' ViC\_ISFC\_ave (62×62)');
subplot(224)
    imagesc(bsl_ISFC_ave); colorbar; axis image; title(' bsl\_ISFC\_ave (62×62)');

%% F2F
t_values = zeros(62);

for m=1:62

    for n=1:62

        edge_1 = squeeze(F2F_ISFC_all(m,n,:));
        edge_2 = squeeze(bsl_ISFC_all(m,n,:));

        % 计算 t 统计量
        [~, ~, ~, stats] = ttest(edge_1, edge_2);

        t_values(m,n) = stats.tstat;
    end

end

[F2F_significant_edges_positive, F2F_significant_edges_negative,F2F_significant_cluster_matrices,F2F_max_pos,F2F_pos_idx,F2F_max_neg,F2F_neg_idx] = scout_significant_edges_test(t_values,F2F_ISFC_all, bsl_ISFC_all);

%% TWE
t_values = zeros(62);

for m=1:62

    for n=1:62

        edge_1 = squeeze(TWE_ISFC_all(m,n,:));
        edge_2 = squeeze(bsl_ISFC_all(m,n,:));

        % 计算 t 统计量
        [~, ~, ~, stats] = ttest(edge_1, edge_2);

        t_values(m,n) = stats.tstat;
    end

end

[TWE_significant_edges_positive, TWE_significant_edges_negative,TWE_significant_cluster_matrices,TWE_max_pos,TWE_pos_idx,TWE_max_neg,TWE_neg_idx] = scout_significant_edges_test(t_values,TWE_ISFC_all, bsl_ISFC_all);
%% ViC
t_values = zeros(62);

for m=1:62

    for n=1:62

        edge_1 = squeeze(ViC_ISFC_all(m,n,:));
        edge_2 = squeeze(bsl_ISFC_all(m,n,:));

        % 计算 t 统计量
        [~, ~, ~, stats] = ttest(edge_1, edge_2);

        t_values(m,n) = stats.tstat;
    end

end

[ViC_significant_edges_positive, ViC_significant_edges_negative,ViC_significant_cluster_matrices,ViC_max_pos,ViC_pos_idx,ViC_max_neg,ViC_neg_idx] = scout_significant_edges_test(t_values,ViC_ISFC_all, bsl_ISFC_all);