in_dir1    = 'D:\桌面\Matlab Working path\truth\F2F'; % F2Fsub*.mat 文件夹
in_dir2    = 'D:\桌面\Matlab Working path\truth\TWE'; % TWEsub*.mat 文件夹
in_dir3    = 'D:\桌面\Matlab Working path\truth\ViC'; % ViCsub*.mat 文件夹
sub_ids   = setdiff(1:30);   
n_nodes   = 62;
n_rand    = 100;
min_w     = 1e-6;
want_kbar = true;               % 是否记录平均度


% ===== 运行并只留下两个结果变量 =====
[F2F_net_all, F2F_net_region] = run_pipeline_with_masks(in_dir1, sub_ids, n_nodes, n_rand, min_w, want_kbar);
[TWE_net_all, TWE_net_region] = run_pipeline_with_masks(in_dir2, sub_ids, n_nodes, n_rand, min_w, want_kbar);
[ViC_net_all, ViC_net_region] = run_pipeline_with_masks(in_dir3, sub_ids, n_nodes, n_rand, min_w, want_kbar);
clear want_kbar sub_ids n_rand n_nodes min_w in_dir
%% =========================== 局部函数（不污染工作区） ===========================
function [net_all, net_region] = run_pipeline_with_masks(in_dir, sub_ids, n_nodes, n_rand, min_w, want_kbar)
    n = numel(sub_ids);

    % ---- 汇总结果容器（全局 + 节点 + 掩膜）----
    net_all = struct();
    net_all.sub_ids      = sub_ids(:);
    net_all.Cw           = nan(n,1);
    net_all.Tw           = nan(n,1);
    net_all.Lw           = nan(n,1);
    net_all.Eglob_w      = nan(n,1);
    net_all.Eloc_w       = nan(n,1);
    net_all.Crand        = nan(n,1);
    net_all.Lrand        = nan(n,1);
    net_all.sigma        = nan(n,1);
    net_all.Qw           = nan(n,1);
    net_all.rw           = nan(n,1);
    net_all.mask_density = nan(n,1);
    if want_kbar, net_all.kbar = nan(n,1); end
    net_all.community    = cell(n,1);
    net_all.mask         = cell(n,1);   % << 新增：每被试的掩膜图（62x62 logical）
    net_all.nodes        = repmat(struct( ...
        'degree',[], 'strength',[], 'Cw_node',[], 'Eloc_node',[], ...
        'BC_node',[], 'EC_node',[], 'P_node',[], 'Z_node',[]), 1, n);

    % ---- 主循环：逐被试 ----
    for ii = 1:n
        sid = sub_ids(ii);
        src = fullfile(in_dir, sprintf('sub%d.mat', sid));
        assert(exist(src,'file')==2, '未找到文件：%s', src);

        % 读原图
        rawW = read_62x62(src, n_nodes);

        % ===== 掩膜 & 度（都基于“原图”）=====
        A_mask = rawW > 0;                 % 掩膜：显著边为1，其余为0
        A_mask(1:n_nodes+1:end) = 0;       % 去自环
        net_all.mask{ii} = A_mask;         % 保存掩膜（logical 62x62）

        k_i = degrees_und(double(A_mask)); % 二值度
        if want_kbar, net_all.kbar(ii) = mean(k_i); end

        % ===== 其余指标：用“清洗后的加权图” =====
        W = sanitize_weight_matrix(rawW, min_w);

        % 加权聚类/传递性
        C_nodes = clustering_coef_wu(W);
        Cw = mean(C_nodes,'omitnan');
        Tw = try_transitivity(W);

        % 路径与效率
        Llen = weight_conversion(W,'lengths');
        D    = distance_wei(Llen);
        [Lw, ~] = charpath(D,0,0);
        Eglob_w = efficiency_wei(W);

        % 局部效率（平均）
        Eloc_nodes = efficiency_wei(W, 2);
        Eloc_w     = mean(Eloc_nodes,'omitnan');

        % 小世界基线（固定掩膜，权值置换）
        [Crand, Lrand] = rand_baseline_same_mask(W, n_rand);
        sigma = (Cw/Crand) / (Lw/Lrand);

        % 模块度 + 社区
        [Ci, Qw] = try_louvain(W);

        % 同配性（强度）
        rw = assortativity_wei(W, 1);

        % 节点级（加权）
        s_i   = strengths_und(W);
        BC_i  = betweenness_wei(Llen);
        EC_i  = eigenvector_centrality_und(W);
        P_i   = participation_coef(W, Ci);
        Z_i   = module_degree_zscore(W, Ci);

        % 写入汇总
        net_all.Cw(ii)           = Cw;
        net_all.Tw(ii)           = Tw;
        net_all.Lw(ii)           = Lw;
        net_all.Eglob_w(ii)      = Eglob_w;
        net_all.Eloc_w(ii)       = Eloc_w;
        net_all.Crand(ii)        = Crand;
        net_all.Lrand(ii)        = Lrand;
        net_all.sigma(ii)        = sigma;
        net_all.Qw(ii)           = Qw;
        net_all.rw(ii)           = rw;
        net_all.community{ii}    = Ci(:);
        net_all.mask_density(ii) = nnz(W)/2 / (n_nodes*(n_nodes-1)/2);

        net_all.nodes(ii).degree    = k_i(:);       % （原图）节点度
        net_all.nodes(ii).strength  = s_i(:);
        net_all.nodes(ii).Cw_node   = C_nodes(:);
        net_all.nodes(ii).Eloc_node = Eloc_nodes(:);
        net_all.nodes(ii).BC_node   = BC_i(:);
        net_all.nodes(ii).EC_node   = EC_i(:);
        net_all.nodes(ii).P_node    = P_i(:);
        net_all.nodes(ii).Z_node    = Z_i(:);
        fprintf('sub%d OK | Cw=%.4f  Lw=%.3f  Eglob=%.4f  Eloc=%.4f  sigma=%.2f  Q=%.3f  r_w=%.3f  k̄=%.1f\n',...
            sid, Cw, Lw, Eglob_w, Eloc_w, sigma, Qw, rw, mean(k_i));
    end

    % ===== 七个脑区索引 & 名称 =====
    idx_regions = { ...
        [21 22 25 26 35 36], ...                            % prefrontal
        [3 4 33 34 37 38 51 52 53 54], ...                  % frontal
        [29 30 41 42 45 46], ...                            % central
        [7 8 9 10 13 14 15 16 27 28 31 32 57 58 61 62], ... % temporal
        [11 12 47 48 55 56 59 60], ...                      % parietal
        [5 6 19 20 23 24 39 40], ...                        % occipital
        [1 2 17 18 43 44 49 50]};                           % limbic
    region_names = {'prefrontal','frontal','central','temporal','parietal','occipital','limbic'};

    % ===== 按脑区把各“节点指标”求平均 → 29×7 =====
    metrics = fieldnames(net_all.nodes(1));
    nR = numel(idx_regions); nS = numel(net_all.sub_ids);
    net_region = struct();
    for m = 1:numel(metrics)
        metric_name = metrics{m};
        mat = nan(nS, nR);
        for si = 1:nS
            vals = net_all.nodes(si).(metric_name);
            for ri = 1:nR
                mat(si,ri) = mean(vals(idx_regions{ri}), 'omitnan');
            end
        end
        net_region.(metric_name) = mat;
    end
    net_region.region_names = region_names;
    net_region.sub_ids      = net_all.sub_ids;
end

% ------------------------- 工具函数（仅在本脚本内可见） -------------------------
function W = read_62x62(matfile, n_nodes)
    S = load(matfile);
    if isfield(S,'PLV') && isnumeric(S.PLV) && isequal(size(S.PLV), [n_nodes n_nodes])
        W = S.PLV; return;
    end
    fns = fieldnames(S);
    for k = 1:numel(fns)
        v = S.(fns{k});
        if isnumeric(v) && isequal(size(v), [n_nodes n_nodes])
            W = v; return;
        end
    end
    error('在 %s 中未找到 %dx%d 数值矩阵。', matfile, n_nodes, n_nodes);
end

function W = sanitize_weight_matrix(W, min_w)
    W(~isfinite(W)) = 0;          % NaN/Inf -> 0
    W(W < 0)        = 0;
    W = (W + W.')/2;
    W(1:size(W,1)+1:end) = 0;
    W(W > 1) = 1;                 % PLV ∈ [0,1]
    small_pos = (W > 0) & (W < min_w);
    W(small_pos) = min_w;         % 防止 1/w 超长路径
end

function Tw = try_transitivity(W)
    try
        Tw = transitivity_wu(W);
    catch
        Tw = NaN;
    end
end

function [Ci, Qw] = try_louvain(W)
    try
        [Ci, Qw] = community_louvain(W);
    catch
        [Ci, Qw] = modularity_und(W);
    end
end

function [Crand, Lrand] = rand_baseline_same_mask(W, n_rand)
    nz = find(W>0);
    Cr = zeros(n_rand,1);
    Lr = zeros(n_rand,1);
    for r = 1:n_rand
        Wr = W;
        Wr(nz) = Wr(nz(randperm(numel(nz))));
        Cr(r)  = mean(clustering_coef_wu(Wr), 'omitnan');
        Llen_r = weight_conversion(Wr,'lengths');
        Dr     = distance_wei(Llen_r);
        [Lr(r), ~] = charpath(Dr,0,0);
    end
    Crand = mean(Cr,'omitnan');
    Lrand = mean(Lr,'omitnan');
end
