# 5 被试教学数据集（demo5）

一套用于学习与测试 SFCX 工具箱的合成数据。它不是真实被试数据，也不用于得出任何科学结论；
它的用途是让第一次使用工具箱的人能在 5 分钟内走通「导入 → 选方法 → 跑分析 → 看结果 → 导出」的完整流程。

所有文件由 `examples/make_example_data.py` 用固定种子 `20261105` 生成，可完全复现：

```bash
python examples/make_example_data.py      # 重新生成
python examples/verify_example_data.py     # 校验并写出 demo5_expected.json
```

---

## 1. 数据规模

| 项目 | 值 |
|---|---|
| 被试数 | 5（`sub-01` … `sub-05`） |
| 脑区数 | 64 |
| 边数 | 2016（上三角） |
| BOLD 时间点 | 120 |
| TR | 2.0 s |
| 网络数 | 6 |
| SC 连接密度 | 24.6 % |

## 2. 文件清单

| 文件 | 内容 | 形状 / 说明 |
|---|---|---|
| `demo5_group.npz` | 全队列打包（推荐，一次导入） | 键：`sc` `fc` `bold` `networks` `hierarchy` `coords` `behaviour` `covariates` `subject_ids` `n_nodes` `tr` |
| `demo5_group.mat` | 同上，MATLAB v7 格式 | 字段与 `npz` 完全一致，MATLAB 直接 `load` |
| `sub-01_sc.csv` … `sub-05_sc.csv` | 单被试结构矩阵 | 64 × 64，对称，对角线为 0 |
| `sub-01_fc.csv` … `sub-05_fc.csv` | 单被试功能矩阵 | 64 × 64，对称，对角线为 1 |
| `sub-01_bold.csv` … `sub-05_bold.csv` | 单被试 BOLD 时间序列 | 64 × 120（脑区 × 时间） |
| `demo5_networks.csv` | 脑区 → 网络标签 | `region,network`，network 取 0–5 |
| `demo5_hierarchy.csv` | 脑区 → 皮层层级 | `region,hierarchy`，1 = 低阶，7 = 高阶 |
| `demo5_coords.csv` | 脑区 → 球面坐标 | `region,x,y,z`，单位球，用于 spin 零模型 |
| `demo5_behaviour.csv` | 被试 → 行为分数 | `subject,score`（示意用，见下方说明） |
| `demo5_covariates.csv` | 被试 → 协变量 | `subject,age_like,motion_like`（示意用） |
| `demo5_expected.json` | 参考结果 | 见第 5 节 |

## 3. 数据结构上的保证

- **BOLD 与 FC 严格自洽**：BOLD 先由目标相关矩阵生成，FC 再由工具箱自带的 `pearson` 管线从 BOLD 算出。
  因此「导入 BOLD」与「导入 FC」得到的是同一组数字（实测最大差值 `0`）。
- **矩阵合法**：SC / FC 均对称、SC 对角线为 0、FC 对角线为 1，SC 密度 24.6 %。
- **网络与层级一致**：`networks` 与 `hierarchy` 逐行对应，层级由网络序号加轻微抖动构成，
  因此变异轴有明确方向，适合演示轴分析。

## 4. 三种导入方式

### 4.1 网页端

1. 启动：`python -m sfcx.webapp --open`
2. 「结构连接 SC」上传 `demo5_group.npz`（服务端会自动取出 `sc`）
3. 「功能连接 FC」上传同一个 `demo5_group.npz`
4. 「网络标签」「皮层层级」「球面坐标」「行为变量」「协变量」同样上传该 `npz`
5. 选择耦合方式（建议先勾 `communicability`、`direct`、`reliance`），估计量全选
6. 置换次数填 `500`，随机种子 `0`，点「运行」

> 也可以只上传 `bold` 那一路：把 `demo5_group.npz` 放到「BOLD 时间序列」槽位，FC 会由 BOLD 现场计算。
> 单被试 CSV 亦可逐槽上传（此时只有 1 名被试，变异轴没有意义，只用于验证格式）。

### 4.2 Python

```python
from sfcx import io, analyze

coh, report = io.load_cohort("examples/data/demo5_group.npz")
print(coh)                       # Cohort(n=5, n_nodes=64, n_edges=2016, ...)

b = analyze(coh, ["communicability", "direct", "reliance"],
            estimators=("regional", "node", "edge"),
            n_perm=500, seed=0)
print(b["sweep"]["r"])
```

### 4.3 MATLAB

```matlab
addpath('matlab');
S = load('examples/data/demo5_group.mat');
cfg = struct('sc', S.sc, 'fc', S.fc, 'networks', S.networks, ...
             'hierarchy', S.hierarchy, 'coords', S.coords, ...
             'behaviour', S.behaviour, 'covariates', S.covariates);
coh = sfcx.buildCohort(cfg);
b = sfcx.analyze(coh, {'communicability', 'direct', 'reliance'}, ...
                 'Estimators', {'regional', 'node', 'edge'}, ...
                 'NPerm', 500, 'Seed', 0);
```

图形界面：`sfcx.app` → 「载入 SC / FC」选 `demo5_group.mat` → 「载入注释」再选同一个文件（会自动识别
`networks / hierarchy / coords / behaviour / covariates`）。

## 5. 参考结果（用于核对导入是否正确）

设置：`n_perm = 500`，`seed = 0`，零模型为 `spin`（因为提供了球面坐标）。

| 耦合方式 | 估计量 | 轴相关 r | 置换 p |
|---|---|---|---|
| communicability | regional | **+0.503** | 0.002 |
| communicability | node | +0.798 | 0.002 |
| communicability | edge | +0.438 | 0.002 |
| direct | regional | +0.335 | 0.050 |
| direct | node | +0.817 | 0.002 |
| direct | edge | +0.322 | 0.004 |
| reliance | regional | +0.494 | 0.002 |
| reliance | node | +0.482 | 0.022 |
| reliance | edge | +0.009 | 0.461 |

只要你的环境与设置相同，数字应当逐位一致（浮点末位可能有 1e-3 级别差异）。

**这套数字本身就是教学点**：同一批数据、同一个轴，换耦合方式或换估计量，结论的强度甚至方向都会变。
`reliance` 的区域级 r 是 +0.49，但边级 r 是 +0.01 —— 这不是程序出错，而是「估计量决定你读到什么」。

## 6. 数据的构造逻辑（便于讲解）

1. 6 个网络中心随机落在球面上，脑区围绕各自中心分布，得到一个层级梯度（网络 0 → 网络 5 由低到高）。
2. SC = 距离核 × 同网络加权，保留每行最近 14 个邻居（约 25 % 密度），再按被试做小幅扰动。
3. 对群体 SC 做扩散核 `exp(-tau·L)`，得到「结构所能预测的功能模式」`C0`。
4. 每个被试的 FC 是 `C0` 与一个被试特异随机相关矩阵的逐区混合，混合权重随层级升高。
   混合式 `diag(a)·C0·diag(a) + diag(b)·R·diag(b)` 仍是合法相关矩阵，所以数据本身没问题；
   而高阶脑区承担更多个体差异，正是变异轴为正的来源。
5. BOLD 用精确相关生成（保留 AR(1) 频谱），因此从 BOLD 反算的 FC 与存储的 FC 完全一致。

## 7. 使用注意

- **行为分数是示意值**：`demo5_behaviour.csv` 中的分数按「个体级耦合强度 + 噪声」构造，
  目的是让行为关联面板有东西可看，不代表任何真实心理测量。5 名被试也无法支撑任何行为推断。
- **p 值只用于演示流程**：5 名被试的置换检验能跑通、能给出 p，但不足以支持统计结论。
- **不要把这套数据与真实分析结果混用**：它只存在于 `examples/` 下，不参与任何项目分析。
