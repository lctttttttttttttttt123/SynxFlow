# T6_design.md — 磷 (P) 吸附-输运模块 · 接口设计蓝图

> **状态**: 接口设计 (2026-06-15)。本文件是开 T6 实现的权威蓝图。**未进实现, 出蓝图后停, 等用户审过才开实现第一步** (同 T5 纪律: 先定接口再填物理)。
> **原则**: 沿用 T5 可插拔框架纪律 —— 公式/参数/速率常数皆配置驱动, kernel 不硬编码任何磷物理数值; 加全新闭合=加 enum case 重编, 换已实现闭合+改参=纯改配置不重编。
> **载体**: 复用 T5 已成的多组 hC 平流框架 (路线甲-正 `_h_flux` 迎风) + bed_k 存量模式 + closures.h 枚举派发风格。
> **参数底稿**: `handoff/T6_literature.md` (长江流域借鉴文献, 无鄱阳本地实测, 正式值留率定)。

---

## 0. 三个已定 fork (用户裁定 2026-06-15)

1. **吸附动力学 = 动力学默认 + 平衡可切换** (enum 派发, closures.h 风格)
   - 默认: **动力学吸附** (准二级 / 一阶速率, 吸附-解吸有速率常数 k_ads, 从配置读)。
   - 可切换: **瞬时平衡** (每步按 Langmuir 等温线重分配两相)。
   - 理由: 文献吸附平衡时间 6~10 h, 与本模型水流时间尺度同量级, 瞬时平衡可能不成立; 默认走动力学更稳, 平衡留 enum 接口。
2. **两相结构 (都要, 非二选一)**
   - **溶解相磷**: 一个独立的**反应性**水相标量场 (带源汇; 区别于 T2 纯被动 (c) 路 —— 那条只守恒不反应)。
   - **颗粒相磷**: 每个泥沙组带一个**吸附磷存量**, 随该组平流/沉降/再悬浮迁移。
   - 两相经 Langmuir 吸附/解吸闭合交换; 颗粒相存量模式**类比 T5 的 bed_k 存量**。
3. **配置化**: Qmax / K / EPC0 / 吸附速率常数 k_ads **全部配置读、每组独立、吸附闭合可替换**; 沿用 T5 的 POD + enum 派发 + 文件缺失向后兼容 (`phosphorus_on=0` → 字节级回退 T5)。

---

## 1. 固定骨架 (要做对)

- **溶解相**: 守恒量 `hPd = h·Pd` (Pd=溶解磷浓度 [mg/L]); 平流复用路线甲-正 (`_h_flux` 迎风, 逐位不动 h/hU), 平流之上挂吸附源汇。
- **颗粒相**: 每泥沙组 k 一个**悬浮吸附磷场** `hPp_k = h·(C_sed_k · q_k)` (q_k=该组单位泥沙吸附磷载量 [mg-P/g]); 随组 k 同一 `_h_flux` 平流; 沉降/再悬浮时与床面吸附磷存量 `bedPp_k` 转移 (类比 hC_k↔bed_k)。
- **两相交换**: 吸附/解吸闭合在每格、每组、每步计算 Pd↔q_k 的净转移, 守恒 (见 §7)。
- **质量守恒判据**: `Σ(hPd·A) + Σ_k(hPp_k·A) + Σ_k(bedPp_k·A) = 常数` (闭箱无净源)。
- **组数 N、磷开关 = 运行时参数**; 颗粒相绑定 T5 泥沙组 (磷依赖 n_sed>0, 见 §8 scope)。

## 2. 可替换接口 (留口子, 不焊死)

- **吸附闭合**: `sorption_mode_id ∈ {0 动力学(默认), 1 瞬时平衡}`; 每个泥沙组独立参数, 从配置读。**先实现动力学(准二级)+平衡两支**; 其余 (Freundlich / 双点位) 预留 enum 位 + TODO。
- **Langmuir 等温线**: `q*(Pd) = Qmax·K·Pd / (1 + K·Pd)` (平衡载量); 物理共识: 细颗粒 Qmax/K 显著高于粗颗粒 (参数体现, 非硬编码)。
- **EPC0 (零净吸附平衡浓度)**: 每组配置量, 决定某格泥沙是磷的源/汇 —— `Pd > EPC0` → 净吸附 (汇), `Pd < EPC0` → 净解吸 (源)。闭合内作平衡参考浓度。
- **Qmax / K / EPC0 / k_ads / q0 (初始载量)**: **全部配置读, kernel 零写死数值** (沿用 T5 §2)。

## 3. 两相场数据结构 + 挂载点

挂在 T5 sediment 框架之上, 全部为**新增场**, 不改 T5 既有场的语义:

| 场 | 类型 | 类比 T5 | 说明 |
|---|---|---|---|
| `hPd` (+ `Pd` remap) | 单个 on_cell Scalar 场 | hC (被动标量) 但**带反应** | 溶解相磷守恒量; 平流 + 吸附源汇 |
| `hPp_k` (k=0..n_sed−1) | per-group on_cell Scalar (device `Scalar**`) | hC_k (悬浮泥沙) | 悬浮颗粒吸附磷; 随组 k 平流/沉降 |
| `bedPp_k` | per-group on_cell Scalar (device `Scalar**`) | **bed_k** | 床面颗粒吸附磷存量; 与 hPp_k 经 E−D 转移 |
| `q_k` (诊断, 可由 hPp_k/hC_k 导出) | — | — | 单位泥沙载量 [mg-P/g], 吸附闭合的状态量; 实现可存场或现算 |

- **POD 参数** `PhosphorusParams` (每组, 设备可用): `Qmax, K, EPC0, k_ads, q0` + `sorption_mode_id` (全局或每组, 见 §4)。
- **全局** `PhosphorusConfig`: `phosphorus_on (默认0)`, `sorption_mode_id`, `vector<PhosphorusParams> groups` (索引对齐 sediment 组), 溶解相初值/边界句柄。
- **挂载**: 颗粒相 P 参数索引严格对齐 sediment 组 (组 k 的 P 吸附在组 k 的泥沙上); n_phos 组数 == n_sed (或子集, 设计取 ==n_sed, 未指定组退化为不吸附 Qmax=0)。

## 4. 吸附闭合抽象 (新 `cuda_phosphorus_closures.h`, 镜像 cuda_sediment_closures.h 风格)

`__device__ inline`, 枚举派发, 零硬编码物理量:
- `langmuir_qstar(const PhosphorusParams& p, Scalar Pd)` → `Qmax·K·Pd/(1+K·Pd)` (平衡载量 helper)。
- `phos_sorption_rate(const PhosphorusParams& p, int mode_id, Scalar Pd, Scalar q_k, Scalar dt)` → 本步 Pd↔q_k 净转移量:
  - `mode_id=0` 动力学 (默认): 准二级 `dq/dt = k_ads·(q*−q)²·sign` 或一阶 `dq/dt = k_ads·(q*−q)` (实现期定一式, 另一式留 enum 子位); q* 由 langmuir_qstar + EPC0 参考给出。**速率常数 k_ads 从配置读**。
  - `mode_id=1` 瞬时平衡: 每步把 (Pd, q_k) 沿 Langmuir 等温线 + 总磷守恒约束**解析重分配**到平衡点 (一步到 q*)。
  - 返回值符号: >0 净吸附 (Pd→颗粒, Pd↓), <0 净解吸 (颗粒→Pd, Pd↑); 死区 (Pd≈EPC0) → 0。
- **正定/守恒**: 转移量 R5 风格截断 (`Pd≥0`、`q_k≥0`、不超 Qmax 饱和); 配套保证 Pd 减少量 == Σ_k 颗粒增加量 (两相严格守恒)。

## 5. 溶解相 + 颗粒相磷的输运 (h/hU 字节级不变, 红线照旧)

复用 T5 路线甲-正多组融合对流, **磷场作纯新增 append 语句骑乘同一 `_h_flux`**:
- 在 WithSediment kernel 内层, 对 `hPd` 与各 `hPp_k` 用**同一 `_h_flux` 迎风**累加 (各自 upwind 浓度), 与 hC_k 同形; 全部 `// [phos]` 标注, **不与 h/hU 算术共享临时量、不交错** → nvcc 不重排 h/hU 的 FMA → **h/hU 字节级不变** (沿用 T5 §7.1 论证)。
- **吸附源汇 = 独立后处理算子** (类比 E−D 与对流分离): 对流+Euler 推进各 P 场后, 调 `cuPhosphorusSorption` (Pd↔hPp_k 交换) + 颗粒相床面转移 (hPp_k↔bedPp_k, 随 sediment E−D 的同比例质量迁移)。
- **算子 (设计命名, 未实现)** 新 `cuda_phosphorus.cu`:
  - `cuPhosphorusSorption(h, hPd, hPps_dev, hCs_dev, pparams_dev, mode_id, n_phos, dt)`: 每格每组经闭合算 Pd↔q_k 转移, 守恒 + R5 截断。
  - 颗粒相随泥沙 E−D 迁移: sediment E−D 已算每组 `dmass`(泥沙质量迁移); P 床面转移按**同一质量分数**搬 hPp_k↔bedPp_k (沉积泥沙携吸附磷入床, 冲刷带回悬浮)。**为不改已验证的 sediment E−D kernel 算术**, 由 sediment E−D 额外输出每组分数迁移场 (additive 输出, 不动 h/hU), P 算子消费之 —— 接口在 §7。
- **主循环插入点** (cuda_flood_solvers.cu, 仅 `phosphorus_on` 内): 每组 P 边界源 → WithSediment 对流 (P 场骑乘) → Euler → τ_b + sediment E−D (输出 dmass 分数) → **cuPhosphorusSorption + 颗粒相床面转移** → Pd/q_k remap (干格清零) → raster。
- **孪生 kernel 同步义务扩展**: T5 §7.1 同步清单新增 `// [phos]` 行; 改对流/Riemann 仍须两 kernel 同源, P append 不参与 h/hU。

## 6. 配置层扩展 (向后兼容)

**推荐: 独立 `input/phosphorus_setup.dat`** (磷完全可选、与 T5 解耦; 文件缺失 → `phosphorus_on=0` → 字节级回退 T5), 分段 key=value, 镜像 `read_sediment_setup()`:
```
$phosphorus_on
1
$sorption_mode      # 0=kinetics(默认) 1=equilibrium
0
$pgroup 0           # 索引对齐 sediment_setup.dat 的 $group 0
Qmax    0.45        # mg-P/g
K       0.7         # L/mg
EPC0    0.02        # mg/L
k_ads   1.0e-4      # 吸附速率常数 (mode=0 用)
q0      0.0         # 初始单位泥沙载量 mg-P/g
$pgroup 1
...
$dissolved          # 溶解相全局
Pd_init  0.0        # 初始溶解磷 (或走 grid-file Pd)
```
- **C++** `read_phosphorus_setup()` → `PhosphorusConfig` (镜像 `read_sediment_setup`: 未知键忽略向后兼容、缺文件返回 phosphorus_on=0)。
- **Python** `InputModel.set_phosphorus_groups([dict,...], sorption_mode=…)` → 写本文件 + 溶解相初值/边界 (复用 IO 标量通道, 不造新格式)。
- **数值=文献占位/示例输入** (见 T6_literature.md), 非建模决策; 正式 Qmax/K/EPC0/k_ads + 鄱阳级配验证阶段填。
- 备选: 并入 sediment_setup.dat 的 $group 段加 P 键 (利用未知键忽略); 取舍记录: 独立文件**胜**在"磷彻底可选、T5 配置零改、缺文件即关", 故推荐独立文件。

## 7. sediment↔phosphorus 接口契约 (E−D 质量分数传递)

- sediment E−D (`cuSedimentErosionDeposition`) 现算每组 `dmass = ed·dt` (泥沙 hC_k↔bed_k 迁移)。T6 需要的是**分数** `f_k = dmass / hC_k` (或迁移绝对量), 供 P 按比例搬运。
- **设计**: 给 sediment E−D 加一个**可选 additive 输出**每组迁移量场 (`phosphorus_on` 时分配, 否则 nullptr); **不改其 h/hU 无关的既有算术**, 仅多写一个输出 → sediment 单独跑 (n_phos=0) 行为字节级不变 (须回归验证)。
- P 床面转移算子读该输出, 按比例移 hPp_k↔bedPp_k; 总磷守恒 = sediment 守恒 + P 比例搬运守恒。

## 8. scope / 非目标

- **磷依赖 n_sed>0** (颗粒相需泥沙载体)。纯溶解-守恒磷 (无泥沙、无反应) 已可用 T2 (c) 路表达, 不在 T6。
- 默认实现: 动力学 (准二级或一阶, 实现期定) + 瞬时平衡两支; Langmuir 等温线; Freundlich/双点位/非黏性磷预留 enum 不实现。
- 沉速: 颗粒相随泥沙组沉速 (T5 张瑞瑾), 不另立; ≤1μm 档絮凝问题继承 T5 (论文交代单颗粒 vs 絮凝)。
- **T6 绝不动 z** (继承 T5 §4): 形态反馈仍 flip-a-switch 扩展点。

## 9. T6 出口闸设计 (沿用 T5 思路: §自洽 + 解析闸, 全过才 close)

**§自洽 (复用 T5 §F 结构)**:
1. 跑通: n_phos 组在合成/case 跑完无 NaN; sm_120 回归。
2. **可插拔**: 仅改 phosphorus_setup.dat (换 sorption_mode / Qmax/K/EPC0/k_ads / n_pgroups) **不重编**, 行为随之变 (贴前后对比)。
3. **守恒**: 闭箱无净源, `Σ(hPd)+Σ_k(hPp_k)+Σ_k(bedPp_k)` 守恒 <0.1% (全精度 backup)。
4. **无回归**: 泥沙各闸仍过; **磷开启后 h/hU 字节级 diff=0** (z 未动); sediment-only (n_phos=0) 对 T5 .so 行为字节级。
5. **n_phos=0 字节级**: 磷关 → 整体对 T5 shipping .so cmp diff=0 (IDENTICAL); WithSediment kernel 寄存器足迹不因 P append 失控污染 (cuobjdump 复核)。

**解析闸 (钉物理正确性, 至少一条; §自洽过不了的)**:
6. **静水吸附平衡闸 (必做)**: 单/多组, 静水 (hU=0), 恒定泥沙 (无沉降/冲刷, τ_b 死区), 给定初始 Pd + q0。
   解析: 吸附达平衡时两相分配对 **Langmuir 等温线** `q*=Qmax·K·Pd_eq/(1+K·Pd_eq)` + 总磷守恒联立解 (Pd_eq 唯一)。判据: 平衡态两相分配相对误差 <1%。验**等温线系数 + EPC0 + 总磷守恒**正确。合成数据, 不进库。
7. **纯解吸一阶速率闸 (可选, 成本低则做)**: 静水, 恒定泥沙带初始载量 q0, Pd 初值 < EPC0 → 解吸; 验 Pd(t) 向平衡的**一阶/准二级速率**对解析衰减/增长曲线 (镜像 T5 沉降柱闸), 抓 k_ads 系数 + 时间积分错误。

## 10. 红线 / 非目标

- **单 GPU only** (多 GPU single_run 不接磷, 留 TODO); **T6 绝不动 master**; **h/hU 数值路径字节级不动** (磷开 == 磷关 == 纯 flood, 硬闸验)。
- **磷参数全配置化, kernel 零硬编码物理量** (Qmax/K/EPC0/k_ads/q0 皆配置读; 仅普适常数/公式身份在码内)。
- 发布边界钩子常开; 真实流域数据 + 级配留仓库外; 合成验证算例入 validation/ 或仓外。
- 物理数值占位/文献缺省 (T6_literature.md) 把框架跑通; 真实标定 (Qmax/K/EPC0 率定 + 鄱阳级配) 留验证阶段。

## 11. commit 粒度 (走钩子, 审过后实现期用)

接口蓝图 (本文件) / 配置读+Python写 (phosphorus_setup.dat) / 两相场+多组融合对流 P 骑乘 / 吸附闭合 (动力学+平衡 enum) + Pd↔hPp_k 交换 + 颗粒相床面转移 + sediment E−D 分数接口 / 解析平衡闸 + 守恒+字节级测试 各一。

---

> **本文件仅蓝图。未写任何实现代码。出蓝图后停, 等用户审过才开 T6 实现第一步。** 红线: 本轮不碰求解器代码 (除本 .md), 不动 master/多 GPU, 不加磷字段到实现, 流场字节级红线照旧。
