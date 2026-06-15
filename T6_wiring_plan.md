# T6_wiring_plan.md — 磷模块接线计划 (路线依据, 等用户裁两处)

> **性质**: 仓库内白名单 .md (同 T2_wiring_plan.md)。**路线决策文档**, 给用户据此裁 §7 路线 + §4 EPC0 口径。
> **状态 (2026-06-15)**: step 1 配置层**已落地** (commit `53c67dd`/`088b036`)。**用户已裁定 (2026-06-15): (A) §7 = 路线① · (B) §4 = 口径 I**, §4/§5/§7 放行, 按 AUTOPILOT 自主推进。裁定细则见文末"决策块 (已裁)"。本文件原为路线决策依据, 现转为实现期权威接线约定。
> **设计权威**: `T6_design.md` (含第1轮审核修订①~④ + 本轮决策侧重审修订 A~D)。**红线照旧**: 单 GPU、不动 master、h/hU 字节级、不动 z、磷参数全配置化、钩子常开。

---

## A. §7 sediment↔phosphorus 接口 —— dmass 传递两路线 (请裁)

**P 颗粒相床面转移要的量**: sediment E−D 每组每格的迁移质量 `dmass`(沉积 <0 / 冲刷 >0), 让吸附磷按比例随泥沙在悬浮↔床面之间搬运 (沉积泥沙携 hPp_k 入床 bedPp_k; 冲刷泥沙带 bedPp_k 回 hPp_k)。

### 已验证 kernel 的事实依据 (cuda_sediment.cu:57–86)

`cuSedimentErosionDepositionKernel` 每格每组 (L66–82):
```
C   = (h>=h_small) ? hc/h : 0          // L69
ws  = sed_settling_velocity(params[k]) // L70   __device__ inline 闭合
ed  = sed_ED_rate(params[k], tb, C, ws)// L71   __device__ inline 闭合
ed  = clamp(ed, -hc/dt, bed/dt)        // L73–76 R5 正定
dmass = ed * dt                        // L77
hCs[k][index]  = hc + dmass            // L78  ← 就地 MUTATE
beds[k][index] = bed - dmass           // L79  ← 就地 MUTATE
```
关键: **dmass 是逐格纯代数率**, 输入全是只读设备场 (`h, tau_b, hc, bed, params, dt`); kernel 唯一副作用是把 `hCs[k]/beds[k]` **就地改写**。→ 重算 dmass 成本≈零 (几个 flop + 同一批已编译闭合), **与 T2 那条昂贵的 `_h_flux` reconstruction+Riemann 完全不同量级** —— 用户的判断成立。

### 路线 ① — P 算子独立重算 dmass, sediment kernel 一字节不动 (**推荐**)

**接法**: 新算子 `cuPhosphorusBedTransfer` 插在 `cuSedimentBedShear` 之后、`cuSedimentErosionDeposition` **之前** (主循环 cuda_flood_solvers.cu, 仅 `phosphorus_on` 内)。它读**与 E−D 完全相同的 E−D 前值** (`h, tau_b, hc_k, bed_k, params_k, dt`)、调**同一批 `sed_settling_velocity`/`sed_ED_rate` 闭合 + 同一 R5 截断**, 重算出**逐位相同的 dmass_k**, 据此移 `hPp_k↔bedPp_k`:
- 沉积 (dmass<0): 悬浮泥沙下沉比例 `f = |dmass|/hc` → `ΔbedPp_k = f·hPp_k`, `hPp_k -= ΔbedPp_k`。
- 冲刷 (dmass>0): 床面泥沙上扬比例 `f = dmass/bed` → `ΔhPp_k = f·bedPp_k`, `bedPp_k -= ΔhPp_k`。

随后 `cuSedimentErosionDeposition` **原封运行**, 把 hc_k/bed_k 按同一 dmass 改写 → 两相同比例迁移、`q_k=hPp_k/hC_k` 不因输运虚变 (接 §C)。

**次序正确性**: 必须跑在 E−D **改写 hc_k/bed_k 之前** (否则读到 post-E−D 的 hc 会重算出不同 dmass)。插入点的 `tau_b`(来自 cuSedimentBedShear) 与 `h`(来自对流) 在 E−D 前后不变 → 重算输入与 E−D 内一致, dmass 逐位相同。

| 维度 | 路线 ① |
|---|---|
| 碰 `cuSedimentErosionDeposition` | **零** (一字节不动) |
| §7 专项硬闸 (碰已验证 kernel) | **不触发** —— 无需 n_phos=0 字节级 / cuobjdump nullptr 复核, **不必为它重编重验 sediment kernel** |
| 成本 | 一次访存型 kernel pass + 逐格几 flop 重算 (近零, 同用户判断) |
| 代价/风险 | dmass+R5 截断**逻辑两处** (E−D 内 + 重算内) → 漂移风险 |
| 漂移缓解 | 物理单源 (两处调**同一** `sed_*` 闭合); 仅 dmass=ed·dt + 3 行 R5 重复; 加**逐格一致性交叉测试** (跑 E−D 抓每格 Δhc; 跑重算; 断言 `dmass_重算 == Δhc == −Δbed` 逐位) 钉死不漂 |

> 注: 若想把 dmass 提成共享 `__device__` helper 消重, **得让 sediment kernel 改走 helper = 碰已验证 kernel**, 反而触发专项硬闸 —— 故路线 ① 保留**有意的、带测试的重复**, 不做共享重构。

### 路线 ② — sediment E−D 加 additive 输出 (退路)

`cuSedimentErosionDepositionKernel` 加可选出参 `Scalar** dmass_out` (phosphorus_on 时分配, 否则 nullptr); 循环内 `if (dmass_out) dmass_out[k][index] = dmass;`。P 算子在 E−D **之后**消费。

| 维度 | 路线 ② |
|---|---|
| dmass 来源 | **单源** (E−D 算一次, 不重复) |
| 碰 `cuSedimentErosionDeposition` | **是** (加 1 出参 + 1 条件写) → 属"碰已验证 kernel" |
| §7 专项硬闸 (T6_design §7 末 / §9.5 提级) | **必触发**: ① n_phos=0 对 T5 .so codegen 字节级; ② cuobjdump 证 nullptr 路径**零 codegen 扰动** (SASS+寄存器足迹不变) |
| 风险 | 加出参进 ABI + predicated store, **极可能改 WithSediment kernel 的 SASS/寄存器分配** (即便运行期 nullptr), 故硬闸 ② 不保证过 → 可能返工或被迫接受 sediment kernel 一个 SASS delta |

### 推荐: **路线 ①**

依据: dmass 是逐格代数率、输入全只读、调同一已编译闭合, 重算近零成本 (≠ T2 的重构+Riemann); 把 P 转移插在**未改动**的 E−D 之前即可逐位复现 dmass —— **从根上消掉"碰已验证 kernel"的风险**, 不触发 §7 专项硬闸、不必重编重验 sediment kernel。唯一代价是 4 行有意重复 (dmass+R5), 由共享闭合 (物理单源) + 逐格 dmass 一致性交叉测试钉死。**STOP: 请裁路线 ① / ②。**

---

## B. §4 吸附闭合 EPC0 ↔ Langmuir q*(Pd) 口径 (请裁)

**收敛成单一平衡参考** (修订 ② / 决策侧重审 B): kernel 内不得有两个能彼此打架的平衡符号源。

### 物理关系 (把 EPC0↔Langmuir 钉死)

默认一阶律 (修订①锁定): `dq/dt = k_ads·(q*(Pd) − q)`, `q*(Pd) = langmuir_qstar(Pd) = Qmax·K·Pd/(1+K·Pd)`。死区 (净转移 0) 落在 `q = q*(Pd)`。

**EPC0 恒等式**: EPC0 = 使当前载量 q 净吸附为零的溶解浓度, 即 `q*(EPC0) = q`。反解 Langmuir:
```
EPC0(q) = q / ( K·(Qmax − q) )          ⟺ 等价地 q*(EPC0)=q
t=0 时: EPC0_0 = q0 / ( K·(Qmax − q0) )
```
→ **EPC0 是载量相关量** `EPC0(q)`, 不是与 q 无关的独立常数。配置里写死的静态 EPC0, 只在"载量恰为 `q*(EPC0_config)`"那一点等于真死区。

### 打架风险 (= 修订① 作废的 (q\*−q)²·sign 坑的根)

若速率同时引用 `q*(Pd)=Langmuir` **和**一个独立静态 EPC0 (如用 `sign(Pd−EPC0)` 定方向), 两平衡源在"载量≠`q*(EPC0)`"处给出**不一致方向** → 死区附近方向错/不连续。这正是平方项抹符号 + 外挂 sign 两源不一致的病根。

### 两口径 (共享不变量: 速率只有一个平衡参考 q\*(Pd); EPC0 与 q0 不能都自由)

- **口径 I — Langmuir 为唯一参考, EPC0 为导出诊断 (推荐)**:
  速率 `= k_ads·(q*(Pd) − q)`, 只引用 `q*(Pd)`。EPC0 不进速率, 仅作**导出诊断** `EPC0(q)=q/(K(Qmax−q))` 报出 (判源/汇用)。
  配置读 `Qmax, K, k_ads, q0`; §6 的 `EPC0` 键降级为**可选一致性校验** (warn if `EPC0_config ≠ q0/(K(Qmax−q0))`, 因 t=0 载量为 q0) 或直接不进速率。
  → 单符号源、零矛盾, 契合 §4 "死区 Pd≈EPC0 → q*≈q → 速率自然趋 0、连续无跳变" (EPC0 涌现而非外加)。
- **口径 II — EPC0 为主测输入, q0 导出**:
  文献多报 EPC0 (非直接 Qmax·K)。取 `EPC0, Qmax, K` 为给定, **导出初始载量** `q0 = q*(EPC0) = Qmax·K·EPC0/(1+K·EPC0)` (令等温线 t=0 自洽过 (EPC0, q0)); 速率**仍只引用 q*(Pd)**。此时 q0 是导出量, §6 的 `q0` 键降级 (或仅作校验)。

> 二者只差**谁导出谁**: 口径 I 由 q0 导 EPC0; 口径 II 由 EPC0 导 q0。**不变量相同**: 速率单参考 `q*(Pd)`, EPC0 与 q0 不并列自由 (避免过定/打架)。

### 推荐: **口径 I** (Langmuir 主参考, EPC0 诊断)

最干净的"单一平衡参考", 直接落实修订②。若用户更想让 EPC0 当配置锚 (文献口径), 取口径 II 亦可 (则 q0 改导出)。
**注**: `T6_literature.md` 现把 EPC0≈0.02 与各组 Qmax/K 与 q0 **并列登记** —— 任一口径下 {EPC0, q0} 必有一者转为导出; 待 §4 放行时按裁定口径把 `T6_literature.md` 的 EPC0↔Langmuir 关系**钉死** (本轮不改 literature, 只在此给结论+依据)。**STOP: 请裁口径 I / II。**

---

## C. 正定性护栏一致性条款 (修订 C; §5 P 平流期落实, 此处单列)

- **hPp_k 必跟 hC_k 走**: 吸附磷不得孤悬在 `hC_k→0` 的格子。任何把 `hC_k` 清零之处 (干格、对流 remap), **同步把 hPp_k 清零** (镜像 T5 R5 / 干湿正定对 hC 的清零)。
- **q_k 除零护栏**: 诊断/重分配里的 `q_k = hPp_k / hC_k` 必加**泥沙趋零护栏** —— `hC_k ≤ ε` (同 h_small 量级干格判据, 镜像 `C=hC/h` 的 R5) 时强制 `q_k=0`、`hPp_k=0`, **绝不做 /hC_k 除法**。
- 与 §3 Qmax=0 惰性组除零安全并列: `q*` 分母 `1+K·Pd≥1` 无奇点; 任何 `/q*` 或 `/hC_k` 归一化短路为零。
- **落实点**: §5 P 平流 (P 场骑乘 `_h_flux` 后) + 每步 `Pd/q_k` remap (干格清零)。本轮不写码, 实现 §5 时落实。

---

## D. 出口闸加严 (修订 D; 测试设计锁定, 此处登记)

§9.7 **纯解吸一阶弛豫闸** 从"可选"提为**必做**, 与 §9.6 静水吸附平衡闸**并列必过**。
解析解 `q(t) = q* + (q0−q*)·exp(−k_ads·t)` (镜像 T5 沉降柱 exp 闸), Pd 钉定隔离速率项。
理由: §自洽五条 (守恒/字节级…) 全过也抓不到 `k_ads` 系数写错 (守恒只管磷没消失、不管吸附速率对不对); 唯这条解析闸钉死 **k_ads 系数 + 量纲 + 时间积分**。本轮不写码, §4/§5 放行后随实现落地。

---

## 决策块 (已裁 — 2026-06-15 用户裁定, §4/§5/§7 放行)

**(A) §7 = 路线① (三不变量焊死; 路线②作废不实现)**:
1. **顺序**: `cuPhosphorusBedTransfer` 插在 `cuSedimentErosionDeposition` **之前**, 读 E−D 同一批 pre-E−D 只读前值 (h, hC_k, bed_k, τ_b); P 算子**只写 hPp_k/bedPp_k, 绝不碰 hC_k/bed_k**。
2. **单源**: dmass 只由 `sed_ED_rate` 闭合 + 同一 R5 算出, 两算子调**同一闭合, 公式零复制** (复用调用不复写公式)。
3. **守闸**: 新增逐格 dmass 一致性交叉测试, 断言 P 算子 dmass == E−D dmass **每格 bit-exact**, 任一格不等即 STOP。
4. **sediment E−D kernel 一字节不动** (不加出参、不改 ABI)。

**(B) §4 = 口径 I**:
1. kernel 速率恒 `k_ads·(q*−q)`, `q* = langmuir_qstar(Qmax, K, Pd)`; **EPC0 不进 kernel 速率** (仅导出诊断 `EPC0(q)=q/(K(Qmax−q))`)。
2. config 主参数 **Qmax / K / q0**; EPC0 可选输入则 `q0 = q*(EPC0) = Qmax·K·EPC0/(1+K·EPC0)` 在 **Python 层**导出; Qmax/K 与 EPC0 同时给且不自洽 → **warning 不静默择一**。
3. `T6_literature.md` 每组同登 {Qmax, K, q0} 与导出 EPC0 (§4 放行后钉)。

**实现范围 (§11 commit 粒度)**: §4 closures → §5 sorption+P 平流 (含 (C) 护栏) → §7 bed-transfer (路线①) → 孪生 kernel 同步 `// [phos]` 行。**出口闸**见 T_progress 顶部 (六闸全过才 close)。
