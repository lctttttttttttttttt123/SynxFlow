# T2_wiring_plan.md — 单 GPU 被动守恒示踪 C 接线方案

> 状态: **路线甲 选定 (2026-06-13)**, 按 AUTOPILOT §A 自主接线。
> 本文件原为选路线用 (§1 代码依据 / §2 两路线对比保留作记录)。下面"决策与执行约束"是落地准则。

## 决策与执行约束 (路线甲)

**选定路线甲**: flood 原 `cuAdvectionMSWEsCartesian` 调用一字不动 (流场字节级不变, 无需重跑 T1d);
新增**只算 hC_advection 的兄弟算子** `cuAdvectionScalarRiderCartesian`, 复用 flood 重构 + 字节相同 Riemann, 仅输出改为 `hC_flux=_h_flux*c_upwind`。

**约束 1 — 复用名副其实 (甲-正: 缓存通量, 非重算)**:
> 修订 (2026-06-14, 用户裁定): 初版兄弟算子**复制重构+重跑 Riemann** 取通量 → case_90 实测 +48% 墙钟, 属**实现走样**
> (开销来自重算, 非"流场零风险")。正解 = 缓存 flood 已算的界面质量通量, 不重算。
- **flood 主 kernel `cuAdvectionMSWEsCartesianKernel` 每步本就算出 `_h_flux` (界面质量通量) 用于更新 h**;
  给它加一个**可选输出参数 `h_flux_cache`** (size 4*ncells, layout `[i*cell_neighbours_length+index]`),
  把 `_h_flux` 落出 (双干面写 0)。**这是"读出中间量", h/hU 的任何算术一字未改** (字节级 diff 验证; 不触红线6)。
- 兄弟 rider `cuTransportScalarRiderCached`: 直接读 `h_flux_cache` 做 `_hC_flux = _h_flux*c_upwind`, **不重跑 reconstruction/Riemann**。
  开销 ≈ 一个 4*ncells 通量数组的显存 + 一次缓存读 + 迎风累加 → 个位数 %。
- flood 原 `cuAdvectionMSWEsCartesian` 7-参版**保留** (多 GPU single_run 用, 内部传 nullptr, 行为不变);
  新增 `cuAdvectionMSWEsCartesianCacheFlux` (传 cache) 供单 GPU run()。
- **同步义务**: rider 与 flood 平流共用同一 `h_flux_cache` → 天然同通量, 无需手工对齐重构。仅依赖主 kernel 的 `[flux-cache]` 落出行。

**约束 1 终态 (c — 2026-06-14, 已落地并通过硬闸)**:
> 缓存版 (commit 3c400c7) 仍 +22.9% (大头是独立 rider stencil pass)。终态把 hC 通量散度**融进主平流 stencil**
> (`cuAdvectionMSWEsCartesianWithTracer` / `cuAdvectionMSWEsCartesianKernel` 加 nullable hC 参数), 复用 kernel 内活着的
> `_h_flux` 就地累加 `_hC_advection`, **消掉独立 rider pass + 缓存数组**。hC 累加为纯新增语句 (新变量 c_this/c_neib/_hC_advection,
> 不与 h/hU 算术交错/不共享临时量), nvcc 未重排 h/hU 的 FMA。
> **硬闸结果**: case_90 h_max/gauge 对纯 flood **字节级 diff=0** (codegen 未漂移); C≡1 max|C-1|=2e-5 (=重算版);
> **运行时 +9.2%** (2450.9s vs 2244s, 单 GPU 36 天); cuobjdump 10x sm_120。并入 (b): 取负折进 Euler(-dt)、删冗余 hC=h*C。
> 缓存版 (3c400c7) 与重算版 (b212942) 留在历史作为 fallback/记录; 终态为融合版。

**约束 2 — 字节级旁证 (比 C≡1 更早抓意外改动)**:
- 接线编译后, 对 `case_90_n020` 跑一发: 因甲不动 h/hU, **h_max/gauge 必须与接线前逐位一致 (diff 全 0 字节级, 不是阈内)**;
  这直接坐实流场未被碰。**同一发**里 C0=1、入流 C=1, 验 `max|C−1|<1e-4` (C≡1 出口闸)。两结果一起进 `T_progress.md`。
- 接线前先存接线前基线产物 (h_max/gauge) 供逐位比对; 该产物在仓库外, 不进库。

**约束 3 — push 节奏**: 每过一个阶段闸 push 一次 (fork 留阶段检查点); push 前自查钩子放行 + 文件清单无数据; 401 → STOP。

**约束 4 — T4 过账后 STOP** 等用户定 T5 方案 (红线3), 不径直进 T5。

---

> (以下为选路线阶段的分析与对比, 保留作记录)
> 原则: **复用输运, 不改水动力**。flood h/hU 流场基线一字不动, 只把 C 作被动标量挂上, 用 flood 每步同一套界面质量通量迎风平流 hC。
> 原则 (用户 2026-06-13 裁定): **复用输运, 不改水动力**。flood 的 h/hU 流场基线 (刚验收的 T1d) 一字不动,
> 只把 C 作被动标量挂上, 用 flood 每步**同一套界面质量通量**迎风平流 hC。**禁止**把平流算子整体换 SRM。

范围: 仅单 GPU `run()`; 多 GPU `single_run()`/`run_mgpus()` 不碰 (留 `// TODO(T2)`)。T2 出口闸 = C≡1 (逐帧 `max|C−1|<1e-4`)。

---

## 1. 复用基线分析 (照搬什么 + 代码依据)

### 1.1 debris 怎么挂 C (参考实现)
`cuda_debris_flow_solver.cu` 的 C/hC 模式 (字段→初始化→主循环→写盘):
| 步 | 行 | 代码 |
|---|---|---|
| 读 C 初值 | 172 | `fvScalarFieldOnCell C_host(.., completeFieldReader("input/field/","C"));` |
| device C | 194 | `cuFvMappedField<Scalar,on_cell> C(C_host, mesh_ptr_dev);` |
| hC | 220 | `cuFvMappedField<Scalar,on_cell> hC(h, partial);` |
| hC_advection | 239 | `cuFvMappedField<Scalar,on_cell> hC_advection(hC, partial);` |
| C_writer | 228 | `cuGaugesWriter<Scalar,on_cell> C_writer(.., C, ".../gauges_pos.dat", "output/C_gauges.dat");` |
| 边界初值 | 336-337 | `C.update_time(t,0); C.update_boundary_values();` |
| **初始化 hC=h·C** | 340 | `fv::cuBinary(h, C, hC, multiply_scalar);` |
| 主循环平流 | 369 | `fv::cuTransportNSWEsSRMCartesian(gravity,h,z,z_gradient,hU, hC, h_advection,hU_advection,hC_advection);` |
| 负号 | 374 | `fv::cuUnaryOn(hC_advection, ... -1.0*a);` |
| **hC Euler 步** | 378 | `fv::cuEulerIntegrator(hC, hC_advection, dt, t);` |
| **改 h 后重映射** | 437-440 | `C=hC/h(divide_scalar); C.update_time/boundary; hC=h·C` |
| **干格清零** | (461 同款) | debris 行 461 `cuBinaryOn(hU,h,filter)`; hC 同款 `b<=1e-10?0:a` |
| C 栅格写盘 | 481 | `raster_writer.write(C, "C", t_out);` |
| C backup | 492 | `cuBackupWriter(C, "C_backup_", t);` |

### 1.2 关键发现 — transport kernel **不是另一套格式**, 是 flood 自身格式 + 被动 hC 骑乘

flood `run()` 平流走 `cuAdvectionMSWEsCartesian` → kernel `cuAdvectionMSWEsCartesianKernel`
(`cuda_advection_NSWEs.cu:836`) → Riemann `cuHLLCRiemannSolverSWEs` (`cuda_riemann_solvers.cuh:65`)。
debris 走 `cuTransportNSWEsSRMCartesian` → kernel (`cuda_transport_NSWEs.cu:118`) → Riemann `cuHLLCRiemannSolverSWEsTransport` (同文件:33)。

**复核结论 (两处 diff):**
1. **`cuHLLCRiemannSolverSWEs` 与 `cuHLLCRiemannSolverSWEsTransport` 字节级相同** (diff 空, 仅函数名不同)。
2. **两 kernel body 仅差"附加 hC 项"** —— diff 全是 `>` 新增行:
   `hC_this`/`c_this`、`hC_neib`/`c_neib` 读入与干格 c=0; `_hC_flux = (_h_flux>=0)? _h_flux*c_this : _h_flux*c_neib`;
   `_hC_advection += _hC_flux*area/volume`; `hC_advection[index]=...`。
   **没有任何一行触及 `_h_flux`/`_hU_flux`/`_z_flux`/`_h_advection`/`_hU_advection`** (那个 flux 调用换成了字节相同的函数)。

→ 即: **h/hU 的数值算法在两 kernel 中完全一致**; transport kernel 内部 `_hC_flux` 复用的正是 `_h_flux` (同一界面质量通量, 1st-order upwind: 出流取本格 c、入流取邻格 c)。这就是用户要的"复用 flood 每步质量通量给 hC"。
唯一非源码层差异: 两个独立编译的函数, FMA 收缩/寄存器调度可能带来 **≤舍入级** 的 h/hU 数值抖动。

### 1.3 干湿正定 (R5, 两条都要)
- transport kernel 内部已处理: 本格/邻格 `h<1e-10 → c=0`, 双干跳过整面。
- 调用方补两道 (debris 同款):
  1. 主循环末 `momentum_filter` 同款对 hC 清零: `cuBinaryOn(hC, h, [](a,b){return b<=1e-10?0:a;})`;
  2. **任何改 h 之处之后立即 `C=hC/h(h>临界否则0); hC=h·C`**。⚠ flood `run()` 的 `cuTotalSourceSink`(行 307)直接改 h,
     必须在其后补此重映射 (被动示踪语义: 降雨=纯水稀释、下渗/抽水=带走纯水, c 不变)。干格 c 强制 0。

### 1.4 Python IO / 边界 — 零新格式, 复用 debris 'C' 通道
recon §2.5/§3.4 已确认: `Boundary.py`/`InputModel.py` 已完整支持 `'C'` 边界源与读写 (debris 在用)。
唯一改动: `InputModel` 的 `__attributes_default` / `__grid_files` / `_file_tag_list` 各加 `'C'` (共 3 行)。**接线阶段再动, 本 STOP 不改。**

---

## 2. 两条候选路线 (请选)

### 路线甲 (优先) — flood 现有框架内加 hC, 完全不动 h/hU 数值路径
**做法**: 保留 flood `run()` 的 `cuAdvectionMSWEsCartesian(...)` 调用**逐字节不变** (h/hU 基线零风险);
新增一个**只算 hC_advection 的兄弟算子** `cuTransportScalarRiderCartesian(gravity,h,z,z_gradient,hU, hC, hC_advection)`
—— 其 kernel = flood kernel 的重构 + 字节相同的 Riemann, 但**只输出 hC_advection**, 不写 h/hU。
hC 的界面通量 `_hC_flux=_h_flux*c_upwind` 用与 flood 完全相同的 `_h_flux` 重算 (同一 reconstruction + 同一 Riemann → 同值)。
- **h/hU**: flood 原调用未改 → **流场字节级不变, 无需重跑 T1d**。
- **代价**: hC 这一步把界面通量再算一遍 (一次额外 advection-kernel 的开销, 被动标量可接受); 新增一个与 flood kernel
  近重复的 hC-only kernel (复用同一算法, 仅改输出选择, 非新算法)。
- **接线点 (run() 内, 单 GPU)**: 字段声明 (§1.1) + 行290后加 `cuTransportScalarRiderCartesian(...)` + 负号 + `cuEulerIntegrator(hC,...)`
  + 改 h 处重映射 (§1.3) + 干格清零 + 写盘 (C_writer/raster/backup) + lib 加新算子 (.cu/.h) + CMake 已含该 src 目录 (无需改)。

### 路线乙 — 退用现成 transport kernel (recon §3.2 原案)
**做法**: 把 flood `run()` 行290 的 `cuAdvectionMSWEsCartesian(...)` **替换**为
`cuTransportNSWEsSRMCartesian(gravity,h,z,z_gradient,hU, hC, h_advection,hU_advection,hC_advection)` (debris 同款), h/hU/hC 一把出。
- **代码最省** (~20 行, 复用现成 kernel, 不新写)。§1.2 证明 h/hU 算法相同, 改动仅舍入级。
- **强制自证闸 (用户要求)**: 接线编译后, **先重跑一次 case_90_n020 完整 36 天, 用 T1d 同口径 reconcile 比 h_max/eta/hU/水量平衡;
  必须满足全部 T1d 阈值 (eta RMSE≤2cm 单站≤5cm、h_max 湿区≤2cm 湿格差≤0.5%、残差率漂移≤0.1pp)。** 不过则 STOP 报你。
- **代价**: h/hU 经一个不同的已编译函数, 理论上有舍入级漂移 (§1.2); 须以上述 T1d 自证背书。

### 2.1 我的判断 (代码依据, 非空结论)
- **甲可行**: 因为 (a) flood kernel 内 `_h_flux` 是界面质量通量, 迎风标量更新只需 `_h_flux*c_upwind` —— transport kernel
  行 88-93 已演示该模式; (b) 复制 flood kernel 的 reconstruction+Riemann 到一个 hC-only 兄弟 kernel 即得同一通量, 算法零改写;
  (c) flood 原 h/hU 调用不动 → **流场基线零风险**, 这是用户"不动水动力"的最强保证。代价是通量重算一遍 + 一个近重复 kernel。
- **甲 vs 乙 的本质**: 二者 h/hU 数学相同 (§1.2 diff 证)。**甲**把 h/hU 留在原函数 (零风险, 多一次通量计算);
  **乙**把 h/hU 也交给 transport 函数 (最省代码, 但须 T1d 自证背书舍入级无害)。
- **倾向**: 若你要"流场绝对一字不动、不想再跑一次 36 天自证" → **选甲**; 若你接受用一次 T1d 自证换取最省代码/零重复 kernel → **选乙**。
  我**推荐甲**: 与你"复用输运不改水动力"的裁定字面一致, 且省掉一次 ~40min 自证跑; 多算一遍界面通量对被动标量是可接受开销。

---

## 3. 两路线共用 (选定后执行, 本 STOP 不做)
- Python: `InputModel` 三处加 `'C'` (§1.4)。
- sm_120: 编译走 `build_env.sh`; build 后 `cuobjdump --list-elf` 验含 sm_120 (回归)。
- commit 粒度: 通读计划(本文件) / 接线 / 编译通过 / C≡1 测试 各一, 走 pre-commit 钩子。
- C≡1 测试在仓库外 (`validation/`), 用 `baseline/case_90_n020_wsl` 流场, 全域 C0=1 + 入流 C=1; 数据不进库。
- 多 GPU: `single_run()` 留 `// TODO(T2): C/hC 多 GPU halo 未接` 注释, 不留半接状态。

## 4. 决策请求
**请选 路线甲 或 路线乙。** 选定后我按 AUTOPILOT 自主推进 (接线→sm_120 编译回归→C≡1 出口闸); 乙则额外先过 T1d 自证闸。
