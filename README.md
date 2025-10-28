# Dyadic Neural Synchronization: Differences  Between Offline and Computer-Assisted Online  Verbal Interaction

本仓库整理了用于 “Dyadic Neural Synchronization: Differences  Between Offline and Computer-Assisted Online  Verbal Interaction” 相关脑电分析的 MATLAB 代码，涵盖从双人被试的跨脑同步（ISC）计算，到显著边的功能网络指标统计。代码主要用于校验与复现论文中的分析流程。

## 目录结构概览

```text
Verbal-Interaction/
├─ 2-6/
│  └─ 04 baseline_and_time_ZPQ.R           # R 语言脚本，计算基线与时间窗口统计量
├─ Major 1-3/
│  ├─ cal_isc.m                            # 主脚本，批量计算多频段 ISC
│  ├─ cal_source_isc.m                     #  ISC 计算
│  ├─ isceeg.m                             # 来自 Cohen et al. (2016) 的 ISC 核心实现
│  └─ ...
├─ Major 1-6/
│  ├─ netmetrics.m                         # 调用 Brain Connectivity Toolbox 计算网络指标
│  ├─ scout_isfc_*.m                       # 映射显著边并输出网络统计的辅助脚本
│  └─ ...
└─ README.md
```

> **提示**：原始数据文件（`.set`、`.mat`、`.xlsx` 等）未包含在仓库中，脚本中出现的 `D:\` 路径为作者本地环境示例，使用时请根据实际存储位置修改。

## 环境与依赖

- MATLAB R2018a 或更高版本。
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php)（`pop_loadset`, `pop_eegfilt`, `topoplot` 等函数）。
- [Brain Connectivity Toolbox (BCT)](https://www.brain-connectivity-toolbox.net/)（`clustering_coef_wu`, `efficiency_wei`, `participation_coef` 等函数）。
- 可选：R（>= 4.0）以运行 `2-6/04 baseline_and_time_ZPQ.R`。

请确保上述工具箱已正确安装并添加到 MATLAB 路径中，BCT 的函数须与 `Major 1-6` 目录下脚本位于同一 MATLAB 搜索路径内。

## 数据准备

1. **EEG 数据**：`Major 1-3/cal_isc.m` 双人被试 `.set` 文件成对存储于 `savefolder` 指定的子目录中（如 `F2F_WF/sub1A.set`, `F2F_WF/sub1B.set`）。请将实际数据路径替换脚本中的示例路径。
2. **通道位置信息**：`locfile` 指向 32/64 通道 `.loc` 文件，用于拓扑绘图与 ISC 计算。
3. **网络矩阵**：`Major 1-6/netmetrics.m` 读取每位被试的 62×62 功能连接矩阵（字段名称可为 `PLV` 或其他 62×62 数组），文件命名约定为 `sub#.mat`。

## 使用指南

### 1. 批量计算 ISC（`Major 1-3`）

1. 打开 `cal_isc.m`，根据实验条件修改 `savefolder`、数据根目录以及输出路径。
2. 确保 EEGLAB 已加载，并能调用 `pop_loadset` 等函数。
3. 运行脚本后，会对 delta、theta、alpha、beta 等频段执行滤波、配对、`isceeg` 计算，并汇总至 `ISC_*_all` 变量。
4. 脚本会记录长度不匹配的被试对到 `erro` 列表，请检查并处理异常数据。

如需在源空间上计算 ISC，可参考 `cal_source_isc.m`，流程与 scalp 版本类似，但输入为源空间的时序矩阵。

### 2. 显著边网络统计（`Major 1-6`）

1. 将显著边阈值处理后的连接矩阵保存为 `sub#.mat`，并更新 `netmetrics.m` 开头的 `in_dir*` 路径。
2. 确保 BCT 函数可被 MATLAB 检索。
3. 运行脚本后将输出：
   - `net_all`：逐被试的全局/节点网络指标、掩膜和社区划分结果；
   - `net_region`：按七个脑区求平均后的节点指标矩阵。
4. `scout_isfc_ncpt2.m`、`scout_isfc_permutest.m` 等脚本提供对显著边的筛选、置换检验及可视化支持，可按需调整参数。

### 3. R 脚本（`2-6`）

该脚本用于基线与时间窗统计量的整理，输入数据（如 `Re_ISC_end_filled_end.xlsx`）用以监测有无残余效应(Residual Effect)，并根据注释修改读写路径。运行方式：

```r
Rscript 2-6/"04 baseline_and_time_ZPQ.R"
```

## 建议的运行顺序

1. 预处理 EEG 数据并导出 `.set` 文件；
2. 运行 `cal_isc.m` / `cal_source_isc.m` 生成 ISC 指标；
3. 将显著边矩阵输入 `netmetrics.m` 计算网络属性；
4. 使用 R 脚本进一步整合统计结果或绘制图表。
