# T5_framework_plan.md — 可插拔泥沙输运框架 (接口设计)

> 原则 (用户裁定): 通用**可插拔**框架, 非鄱阳湖算例。公式/参数/组数皆配置驱动, kernel 不硬编码任何物理数值。
> 验收: 换河/换参数/换公式 → 只改配置不改 kernel 源码就能跑。鄱阳湖只是验证案例。
> 接口 3 fork 已拍板 (2026-06-14): 枚举派发 / 分段 .dat+Python / vector<field>+融合多组 stencil。
> 物理数值用占位/文献缺省把框架跑通; 真实鄱阳湖参数 (两组 D50、张瑞瑾沉速、τ) 留验证阶段填入配置。

---

## A. 配置格式 — `input/sediment_setup.dat` (分段 .dat, C++ 简单解析, Python 写)

```
$n_groups
2
$group 0
settling_id   1            # 0=const  1=Zhang(张瑞瑾)  2=Stokes  3=floc(絮凝修正)
closure_id    0            # 0=Partheniades 双阈值(黏性)  [1=Shields 非黏性: 预留, 暂不实现]
cohesive      1
rho_s         2650.0
D50           6.0e-6
w_s           5.0e-5       # const 沉速 / 或作 Zhang 的回退值
tau_ce        0.15         # 临界冲刷应力 (Pa)
tau_cd        0.06         # 临界淤积应力 (Pa)
M             2.0e-4       # Partheniades 冲刷系数 (kg/m2/s)
$group 1
settling_id   1
closure_id    0
...
```
- **所有物理量从此文件读, kernel 内零写死数值。** 缺省值=文献占位, 非建模决策。
- Python: `InputModel.set_sediment_groups([{settling_id:.., closure_id:.., rho_s:.., D50:.., w_s:.., tau_ce:.., tau_cd:.., M:.., cohesive:..}, ...])` → 写出本文件 + 每组初始浓度 field `input/field/C0`, `C1`, ...（沿用 'C' grid-file 机制, 扩成 C{k}) + 每组边界 `C{k}_BC_*.dat` (复用 'C' 通道)。
- 求解器侧: 一个 `SedimentConfig` 结构体 vector (每组一条), `read_sediment_setup("input/sediment_setup.dat")` 填充。

## B. 闭合抽象 — 枚举派发 (运行时选, 编译期实现)

```c
// __device__ 内联, switch 派发; 所有参数从 SedimentParams (每组) 传入, 无写死
__device__ Scalar settling_velocity(int settling_id, const SedimentParams& p, Scalar h, Scalar C, ...);
__device__ Scalar ED_rate_closure(int closure_id, const SedimentParams& p, Scalar tau_b, Scalar h, Scalar C, Scalar w_s);
```
- `settling_id`: 0=const(p.w_s)  1=Zhang 张瑞瑾  2=Stokes  3=floc。**先实现 const + Zhang**, 其余 case 留 stub+TODO(返回 const 回退并 printf 一次)。
- `closure_id`: 0=Partheniades 双阈值黏性:
  - 淤积 (τ_b<τ_cd): `D = w_s * C * (1 - τ_b/τ_cd)`;  冲刷 (τ_b>τ_ce): `E = M*(τ_b/τ_ce - 1)`;  中间带=0。
  - `ED_rate = E - D` (>0 净冲刷上浮, <0 净淤积)。
  - `closure_id=1` (Shields 非黏性) 预留枚举位, **本研究不实现** (留扩展点)。
- **选已实现公式 + 改参数 = 纯改配置不重编; 加全新公式 = 加一个 case 重编** (符合用户"可插拔"的现实含义)。

## C. 组数 N + 多组对流 — vector<field> + 融合多组 stencil

- N 运行时: `std::vector<cuFvMappedField<Scalar,on_cell>> C(N), hC(N), hC_adv(N), bed(N)` (bed=床面可冲刷沉积量/组)。
- **对流复用 路线甲-正**, 但泛化为多组: 新增 `cuAdvectionMSWEsCartesianWithSediment(... , Scalar** hC, Scalar** C_bound, Scalar** hC_adv_out, int n_groups)`;
  kernel 在算出 `_h_flux` 后**内层 loop k=0..N-1** 用同一 `_h_flux` 做迎风累加各 `hC_adv[k]` —— 一趟 stencil 出 h/hU + 全部 N 组。
  - hC 累加仍为**纯新增语句**(不与 h/hU 算术交错), 保 h/hU 字节级不变 (硬闸照旧验)。
  - 设备侧 `Scalar** hC` = 主机组装的指针数组 (cudaMemcpy 到 device 一个 Scalar*[N])。
- T2 单标量 passive tracer `WithTracer` 保留不动 (已验证); 多组 sediment 走新 `WithSediment` 路径 (N≥1)。

## D. 与水流耦合 (要对)

- **床面切应力** τ_b 取自 flood Manning 摩阻: `τ_b = ρ_w * g * n² * |u|² / h^(1/3)` (|u|=|hU|/h, 干格 0)。一个 helper kernel 算 τ_b 场。
- 源汇挂对流之上 (每组每步): 算 τ_b → 每组 `ED_rate_k = closure(τ_b, C_k, w_s_k, params_k)` → **R5 正定截断** (`ED_rate ≥ -hC_k/dt` 不抽负; 冲刷 `≤ bed_k*(...)/dt` 不超可冲刷量) → Euler `hC_k += ED_rate_k*dt` → 同步 `bed_k -= ED_rate_k*dt` (质量守恒: 悬浮+床面) → `C_k=hC_k/h` → 干格清零。
- 沉降/再悬浮通过 ED_rate 反馈到浓度 (沉降=负 ED 进 bed, 再悬浮=正 ED 出 bed)。

## E. 不碰 / 红线
- 单 GPU only (多 GPU `single_run` 不接, 留 TODO)。h/hU 数值路径字节级不动 (硬闸验)。发布边界钩子常开。真实流域数据留仓库外。

## F. T5 出口闸 (建议, 待用户确认)
1. **跑通**: N=2 占位组 (文献缺省) 在合成/case_90 流场跑完, 无 NaN, sm_120 回归。
2. **可插拔验证**: 仅改 `sediment_setup.dat` (换 closure_id/settling_id/参数/n_groups), **不重编** kernel, 行为随之变 (贴前后对比)。
3. **守恒**: 闭箱无净源时, 每组 总悬浮+床面质量 Σ(hC_k·A)+Σ(bed_k·A) 守恒 <0.1% (全精度 backup)。
4. **耦合 sanity**: 高切应力格冲刷 (C 升)、低切应力格淤积 (C 降向 bed); 量纲守恒; 平衡态自洽。
5. **无回归**: T2 passive C≡1 仍 <1e-4; 纯 flood h/hU 字节级 diff=0。

## G. commit 粒度
接口 spec(本文件) / 配置读+Python写 / 多组对流 kernel / τ_b+E-D 闭合 / 测试 各一, 走钩子。
