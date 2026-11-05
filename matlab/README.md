# SFCX MATLAB 工具箱

与 Python 端同一套耦合定义、变异估计、零模型与轴统计量；界面与网页版一一对应。

## 快速开始

```matlab
addpath('<repo>\matlab')            % <- path to this matlab folder
sfcx.app                    % 图形界面
sfcx.describe()             % 打印定义 / 估计量 / 零模型清单
```

无界面路线：

```matlab
coh = sfcx.demo(60, 40, 0);                  % 合成队列，便于先跑通
b = sfcx.analyze(coh, {'communicability', 'direct'}, ...
                 'estimators', {'regional', 'node', 'edge'}, ...
                 'nPerm', 500, 'seed', 0);
b.sweep.r        % 定义 × 估计量的轴相关系数
b.metrics(1).estimators(1).axis    % 单格：r / p / n / 零模型
b.tables.axis_r  % CSV 文本
```

## 数据接入

```matlab
coh = sfcx.buildCohort(struct( ...
        'sc', sc, 'fc', fc, 'networks', net, 'hierarchy', hier, ...
        'coords', xyz, 'behaviour', y, 'covariates', cov));
```

| 字段 | 形状 | 说明 |
|---|---|---|
| `sc` / `fc` | `(N, E)`、`(N, P, P)`、`(P, P)` | 上三角展开；矩阵输入自动展开 |
| `bold` | `(N, P, T)` 或 `(P, T)` | 由 `pearson` / `precision` 转成 FC |
| `networks` | `(P,)` | 网络标签 |
| `hierarchy` | `(P,)` | 皮层层级位置 |
| `coords` | `(P, 3)` | 球面坐标，启用 spin 零模型 |
| `behaviour` / `covariates` | `(N,)` / `(N, k)` | 行为与协变量 |

单独调用：

```matlab
fc  = sfcx.boldToFc(bold, 'precision');     % (N, E)
coh = sfcx.buildCohort(struct('bold', bold, 'sc', sc));
```

## 函数对照

| Python | MATLAB | 作用 |
|---|---|---|
| `sfcx.coupling(coh, m)` | `sfcx.coupling(coh, m)` | 一个耦合定义 → 边级或区域级场 |
| `sfcx.regional_readout` | `sfcx.regionalReadout` | 边级场 → 区域级读出 |
| `analysis.regional_variability` | `sfcx.regionalVariability` | 区域级场的跨被试 SD |
| `analysis.variability` | `sfcx.variability` | 边级 / 节点投影的 SD |
| `analysis.axis` | `sfcx.axis` | 与层级轴的 Spearman r 与零模型 p |
| `analysis.layers` | `sfcx.layers` | 网络内 vs 网络间 |
| `analysis.epicentres` | `sfcx.epicentres` | 焦点区域与分离度 z |
| `analysis.deviation_spectrum` | `sfcx.deviationSpectrum` | 个体偏差谱四尺度斜率 |
| `analysis.reliability` | `sfcx.reliability` | split-half 一致性 |
| `engine.analyze` | `sfcx.analyze` | 多定义 × 多估计量批量分析（n≥2） |
| `engine.magnitude` | `sfcx.couplingValues` | 耦合值计算（n≥1，无变异统计） |
| `engine.export_coupling` | `sfcx.exportCouplingValues` | 逐被试 CSV / MAT 导出 |
| `io.batch_stack`* | `sfcx.batchStack` | 文件夹批量堆叠为 N×P×P |
| `custom.register_custom` | `sfcx.registerCustom` | 注册自定义耦合 |
| `io.build_cohort` | `sfcx.buildCohort` | 统一数据接入 |
| — | `sfcx.exportResults` | 导出 CSV / 图 / 报告 |

*Python 端批量导入通过上传多个文件 / npz 完成对应功能。

## 自定义耦合

```matlab
src = sprintf(['function out = couple(sc, fc, iu, iv, nNodes)\n' ...
               '    out = sc .* fc;\n' ...
               'end\n']);
sfcx.registerCustom('my_coupling', src, 'elementwise product');
out = sfcx.coupling(coh, 'my_coupling');

sfcx.listCustom();            % 查看已注册项
sfcx.unregisterCustom('my_coupling');
```

签名 `couple(sc, fc, iu, iv, nNodes)`；返回 `(N, E)`、`(E,)`、`(N, P)`、`(P,)`、`(P, P)` 或标量。也接受函数句柄 `sfcx.registerCustom('name', @myCoupling)`。注册前会在小样本上试算，返回全 NaN 或尺寸不合法会直接拒绝。

## 界面说明（v0.3 重写版）

`sfcx.app` 是四个顶层页签的窗口，中英双语（右下角 EN/中文 切换）：

| 页签 | 内容 | 被试要求 |
|---|---|---|
| **① 数据导入** | 单被试模式（SC + FC/BOLD 逐槽位）或多被试批量模式（选文件夹每人一个文件，或单个含三维数组/元胞的 .mat）；右侧注释槽位（网络、层级、坐标、行为、协变量）与「构建队列」 | 1 名或多名 |
| **② 变异分析** | 定义多选 → 估计量、零模型、置换次数 → 运行；结果子页签：连接轴表、变异剖面、层级与震中、坐标扫描 | **≥2**（<8 会警告噪声大） |
| **③ 耦合值计算** | 定义多选 → 计算每个被试的耦合值；查看组平均矩阵热图与区域值热图；导出 CSV（逐被试边值/区域值）与 MAT（含平均矩阵） | **≥1**（单被试可用） |
| **④ 日志** | 运行记录 + 变异分析结果导出（CSV/MAT） | — |

单被试数据不会报错卡死：页签②会明确提示"变异度统计在数学上未定义，请用页签③"，页签③正常出结果。

### 批量导入的两种方式

1. **文件夹**：每人一个 `.csv/.txt/.tsv/.mat`（SC 批量为 P×P 矩阵；BOLD 批量为 P×T），文件名排序即被试顺序，`sfcx.batchStack` 负责堆叠与尺寸校验；
2. **单个 .mat**：含一个三维数组（N×P×P 或 N×P×T）或元胞数组（每格一名被试）。

### 耦合值计算（命令行）

```matlab
vals = sfcx.couplingValues(coh, {'direct', 'communicability'});
vals(1).meanMatrix      % P×P 组平均矩阵（单被试时即本人矩阵）
vals(1).regional        % N×P 区域级读出
sfcx.exportCouplingValues(vals, outdir, {'csv', 'mat'});
%   <metric>_coupling_edges.csv     逐被试全部边
%   <metric>_coupling_regional.csv  逐被试全部脑区
%   <metric>_coupling.mat           场 + 区域值 + 平均矩阵 + 被试标签
```

区域级定义（Granger、传递熵、RNC/TRR 等）没有 P×P 矩阵，只有区域值——界面与 README 均有说明。

## 数值一致性

```matlab
run('...\matlab\tests\run_tests.m')        % 核心 9/9
run_gui_tests                              % GUI + 批量 + 耦合值导出 6/6
```

覆盖 Cohort、目录、耦合定义、轴与零模型、层级、批量分析、自定义函数、BOLD 转换与 Python 参考值对照，当前 9/9 通过；GUI 测试覆盖 batchStack、三维批量建队、单/多被试耦合值、CSV+MAT 导出、单被试拒绝逻辑与 GUI 构建，当前 6/6 通过。

对照策略（写在 `reference.mat` 里，由 `matlab/tests/make_reference.py` 生成）：

- 12 个快定义 × 3 估计量 = 36 格，其中 **33 格逐格对照**；
- 10 个定义的轴相关系数 |Δr| ≤ 5e-7（机器精度级一致）；
- `reliance` 边级读出有 ~1e-2 量级的浮点累积差异，容差 5e-2；
- `search_info` **排除**：最短路径长度并列时，scipy 与 MATLAB 选的路不同，
  累积搜索信息因此不同——两个运行时内部各自自洽，但不可逐格对照。

参考矩阵改动后需在同一项目环境下重跑 `python matlab/tests/make_reference.py`。
注意：两个运行时的随机流不同，**p 值不应逐一比对**；可比的量是确定性的 r、SD、斜率与 split-half。

## 环境提示

- 组件使用 `uigridlayout` 布局，不依赖 `Units`/`Position`，避免不同 MATLAB 版本的属性差异。
- 若启动时提示某个第三方工具箱与内置函数同名（如 `narginchk`），那是 MATLAB 路径里其他工具箱造成的，与本工具箱无关；注册自定义函数时已临时关闭该警告。
- 大样本建议先跑区域级估计量：边级变异在 360 脑区时是 64,620 维，节点投影与矩阵类定义会更慢。
