# T5_design.md — 可插拔泥沙输运框架 · 执行蓝图

> **状态**: 接口设计四 fork 全定 (2026-06-14)。本文件是开 T5 实现的权威蓝图。**未进实现, 等用户明确"开 T5"。**
> **原则**: 通用**可插拔**框架, 非鄱阳湖算例。公式/参数/组数皆配置驱动, kernel 不硬编码任何物理数值。
> **验收**: 换河/换参数/换公式 → 只改配置不改 kernel 源码就能跑。鄱阳湖只是验证案例之一。
> (超越/取代 T5_framework_plan.md)

## 0. 四个已定 fork
1. **闭合抽象 = 枚举派发**: 每组 settling_id/closure_id (int, 配置读), kernel 内 `switch` 派发到已编译公式; 全参数从配置读。选已实现公式+改参=纯改配置不重编; 加全新公式=加 case 重编。
2. **配置格式 = 分段 .dat + Python 写**: `input/sediment_setup.dat` (§5); `InputModel.set_sediment_groups([...])` 写出。
3. **组数 N + 对流 = vector<field> + 融合多组 stencil**: 运行时 N 个浓度场; 路线甲-正 kernel 泛化为内层 loop k 复用同一 `_h_flux` 一趟出全部组 (§7)。
4. **悬浮输运 only + Δz_accum 形态扩展点 (§4) + 出口闸 = §F 五条 + 解析闸 (§9)**。

## 1. 固定骨架 (要做对)
- 每泥沙组一个浓度场 C_k / 守恒量 hC_k; 对流**复用路线甲-正** (融合 `_h_flux` 迎风, 逐位不动 h/hU)。
- 源汇 E−D 挂对流之上; 干湿正定沿用 R5 (改 h 后 C=hC/h、干格清零; E−D 后 hC 正定截断)。
- 与水流耦合: τ_b 自 flood Manning 摩阻; 冲淤通过 hC_k↔bed_k 转移体现 (§8)。
- **组数 N = 运行时参数, 非编译期常数**。

## 2. 可替换接口 (留口子, 不焊死)
- 沉速: settling_id ∈ {0 const, 1 Zhang(张瑞瑾), 2 Stokes, 3 floc}; 每组独立; 从配置读。**先实现 const+Zhang**, 余 stub+TODO。
- E−D 闭合: closure_id=0 Partheniades 双阈值黏性 (默认实现); closure_id=1 Shields 非黏性 **预留枚举位不实现** (扩展点)。
- τ_ce/τ_cd/M/w_s/D50/ρs/cohesive: **全部配置读, kernel 零写死数值**。

## 3. (并入 §0/§2)

## 4. 形态动力 = flip-a-switch 扩展点 (T5 绝不动 z)
**T5 阶段绝不改 z** —— z 一变 h/hU 变, 路线甲+(c) 的"流场字节级 diff=0"立即作废。框架未验证扎实前不引形态反馈。
- 冲淤通过每组 **bed_k (床面记账存量场)** 与 hC_k 间转移体现; **质量守恒 = Σ(hC_k·A)+Σ(bed_k·A) 闭合, 全程不碰 z**。
- **Δz_accum 计算但不施加**: 每步累计净冲淤的床面高程增量 `Δz_accum += Σ_k (Δbed_k / ρs_k / (1−porosity))`,
  **算出来并存为一个场**, T5 不拿它改 z (computed-but-unapplied)。
- 配置留 `morphology_on` 标志 (默认 false); 注释写明: true 分支 = "把 Δz_accum 加到 DEM z 上并触发水流重算", **未实现的扩展点, 留 TODO**。
- **z 被读进水动力的耦合点**在代码注释明确标注, 使将来接形态动力只需接这一根线 (开关 + 接线, 非动骨架)。

## 5. 配置 `input/sediment_setup.dat` (分段 key=value)
```
$n_groups
2
$morphology_on
0
$group 0
settling_id   1
closure_id    0
cohesive      1
rho_s         2650.0
porosity      0.4
D50           6.0e-6
w_s           5.0e-5
tau_ce        0.15
tau_cd        0.06
M             2.0e-4
$group 1
...
```
- C++ `read_sediment_setup()` → `std::vector<SedimentParams>` (+ morphology_on)。Python `set_sediment_groups([dict,...])` 写本文件 + 每组初始浓度 field `C0`,`C1`,.. (扩 'C' grid-file) + 每组边界 `C{k}_BC_*.dat` (复用 'C' 通道)。
- **数值=文献占位/示例输入, 非建模决策**; 真实鄱阳湖参数 (两组 D50=6μm/1μm、张瑞瑾沉速、1μm 絮凝 vs 单颗粒、τ 来源) 验证阶段填。

## 6. 闭合 (枚举派发, __device__ 内联)
- `settling_velocity(settling_id, p, ...)`: 0→p.w_s; 1→张瑞瑾; 2→Stokes; 3→floc。先 const+Zhang, 余回退 const+一次性 printf。
- `ED_rate_closure(closure_id=0 Partheniades, p, τ_b, C, w_s)`:
  - 淤积 (τ_b<τ_cd): `D = w_s·C·(1−τ_b/τ_cd)`; 冲刷 (τ_b>τ_ce): `E = M·(τ_b/τ_ce−1)`; 中间带 0。
  - `ED_rate = E − D` (>0 净再悬浮, <0 净淤积)。

## 7. 多组融合对流 (h/hU 字节级不变)
- `cuAdvectionMSWEsCartesianWithSediment(..., Scalar** hC, Scalar** C_bound, Scalar** hC_adv_out, int n_groups)`:
  kernel 算出 `_h_flux` 后**内层 loop k** 用同一 `_h_flux` 迎风累加各 `hC_adv[k]`。hC 累加为**纯新增语句** (不与 h/hU 算术交错), 保 h/hU 字节级不变。
- 设备侧 `Scalar**` = 主机组装指针数组拷到 device。T2 单标量 `WithTracer` 保留 (N=1 passive 用)。

### 7.1 KERNEL 同步清单 (twin kernel 维护义务 — 改一个必须同步另一个)
step 2 放行时触发了预案 (b) **分叉**: 单 kernel 泛化让 n_sed=0 也带 MAXSED=8 栈数组 → 寄存器压力
(实测 +26% vs (c) +9.2%), 污染了已验证的被动路径。故保留两个**孪生** kernel:

| kernel | 路径 | REG/STACK (sm_120) | 用途 |
|---|---|---|---|
| `cuAdvectionMSWEsCartesianKernel` (= 路线甲-正 (c)) | n_sed=0 | **REG:64 / STACK:0** | `cuAdvectionMSWEsCartesian` (single_run, nullptr×3) + `cuAdvectionMSWEsCartesianWithTracer` (flood passive 标量 C) |
| `cuAdvectionMSWEsCartesianWithSedimentKernel` | n_sed>0 | REG:72 / STACK:160 | `cuAdvectionMSWEsCartesianWithSediment` (h/hU + 被动 C + N 组泥沙, 一趟) |

**铁律**: 两 kernel 的 **重构 (MSWEs surface reconstruction) + HLLC Riemann + h/hU/_z_flux 算术**必须**逐行同源**。
- 任一方改对流/Riemann 逻辑, **另一方必须同步**, 否则 n_sed=0 与 n_sed>0 的流场会分叉。
- 泥沙相关语句在 WithSediment kernel 里全部是 `// [sed]` / `// [tracer]` 标注的**纯新增 append 语句**
  (新变量、不与 h/hU 算术共享临时量、不交错), 故 nvcc 不重排 h/hU 的 FMA → h/hU 字节级不变 (3c [a] 实证: 泥沙开关 h/hU 逐位 IDENTICAL)。
- **理想**: 将公共重构+Riemann 抽成一个 `__device__` inline 函数两 kernel 共用 (消除复制); 当前为复制+本清单约束 (退一步)。重构若抽公共函数, 本表保留作回归基准。
- 代码内已在两处 (`.cu` kernel 头注释 + `.h` 声明注释) 留了"SYNC: 见 T5_design.md kernel 同步清单"指针。

**验证状态 (step 2, 2026-06-14)**: 3a n_sed=0 对 t2_prewire 字节级 diff=0; 3b cuobjdump 证 n_sed=0 kernel REG:64/STACK:0 无污染;
3c 合成两组对流-only: max|C0−1|=0、max|C1−0.5|=0 (独立)、泥沙开 h/hU 字节级 IDENTICAL。源汇 E−D 仍未加 (step 3)。

## 8. τ_b 耦合 + E−D + R5 + bed 守恒 (每组每步)
1. `τ_b = ρ_w·g·n²·|u|² / h^(1/3)` (|u|=|hU|/h, 干格 0) — helper kernel 出 τ_b 场。
2. 对流+Euler 推进 hC_k (§7)。
3. 每组: `ED_rate_k = closure(τ_b, C_k, w_s_k, p_k)` → **R5 正定截断** (`≥ −hC_k/dt`; 冲刷 `≤ bed_k·(...)/dt` 不超可冲刷) → `hC_k += ED_rate_k·dt`; **同步 `bed_k −= ED_rate_k·dt`** (悬浮+床面守恒); 累计 `Δz_accum` (§4, 不施加)。
4. `C_k=hC_k/h`; 干格 hC_k/C_k 清零 (R5)。

## 9. T5 出口闸 = §F 五条 + 解析闸 (六条全过才 close)
**§F 五条 (自洽)**:
1. 跑通: N=2 占位组在合成/case_90 跑完无 NaN; sm_120 回归。
2. 可插拔验证: 仅改 sediment_setup.dat (换 closure_id/settling_id/参数/n_groups) **不重编**, 行为随之变 (贴前后对比)。
3. 守恒: 闭箱无净源, 每组 Σ(hC_k·A)+Σ(bed_k·A) 守恒 <0.1% (全精度 backup)。
4. 耦合 sanity: 高 τ_b 冲刷 (C↑)、低 τ_b 淤积 (C↓向 bed); 量纲守恒。
5. 无回归: T2 passive C≡1 仍 <1e-4; **泥沙开启后纯 flood h/hU 字节级 diff=0** (z 未动硬证)。

**解析闸 (钉物理正确性, §F 自洽过不了的)**:
6. **纯沉降柱 (必做)**: 单组, 静水 (hU=0), 关冲刷 (E=0, τ_b=0<τ_cd → D=w_s·C), 1D/单柱。
   解析 `C(t)=C0·exp(−w_s·t/h)`。判据: 逐时刻相对误差 <1% (或与一阶时间积分截断误差一致的容差, 取严)。
   验沉降项**系数+量纲+时间积分数值正确**, 不只守恒。合成数据, 不进库。
7. **纯冲刷增长率 (可选, 实现成本低则做)**: 静水, 恒定 τ_b>τ_ce, 初始 C=0; 验 `E=M(τ_b/τ_ce−1)` 线性增长率对解析积分, 抓 E 系数错误。
   成本高则先只做沉降柱, 冲刷靠耦合 sanity + 守恒兜底。

## 10. 红线 / 非目标
- 单 GPU only (多 GPU single_run 留 TODO, 不接 sediment)。**T5 绝不动 z** (§4)。h/hU 数值路径字节级不动 (硬闸验)。
- 发布边界钩子常开; 真实流域数据留仓库外 (合成验证算例入 validation/synthetic/ 或仓库外)。
- 物理数值占位/文献缺省把框架跑通; 真实标定留验证阶段。为 T6 磷模块载体铺垫 (hC_k 多组分 → 磷粒径依赖吸附挂组上)。

## 11. commit 粒度 (走钩子)
接口蓝图(本文件) / 配置读+Python写 / 多组融合对流 kernel / τ_b+E−D 闭合+bed+Δz_accum / 解析+守恒测试 各一。
