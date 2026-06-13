# PUBLICATION_BOUNDARY.md — 这个 fork 能/不能进什么

> **这是公开发布边界的权威说明。** 本 fork 以 **GPL** 衍生方式公开在 GitHub。
> 推到公开 fork 的任何提交都是**不可逆的**: git 历史会永久保留, 即便后续删除/改写,
> 也已被 GitHub 缓存、镜像站、克隆者、网络爬虫抓走。**漏放一次 = 真实流域数据
> public + GPL 永久可见。误拦只浪费十秒。** 故守卫采 **default-deny (白名单优先)**。

强约束执行点: `scripts/pre-commit` 钩子 (版本库内权威副本) + 根 `.gitignore` (软层)。

---

## 1. 放行三条件 (须同时满足)

一个暂存文件**只有**同时满足下面三点才放行, 否则 `pre-commit` 报命中清单并 `exit 1`:

1. **命中白名单** (扩展名 / 具名文件 / 授权路径) —— 见 §2;
2. **不命中关键词黑名单** (鄱阳湖流域专名 + 数据 token) —— 即便扩展名在白名单里也拦, 见 §3;
3. **单文件 ≤ 5 MB**。

附加 (上一轮约定保留): **脚本 (`.sh`/`.py`) 内容不得含驱动器绝对路径** (`F:\` 或 `/mnt/<drive>/`),
防止把"读 F 盘鄱阳湖基线"的脚本提交进来。

> 这不是冗余: 白名单挡住扩展名层 (`.asc/.dat/.nc` 等数据格式), 关键词黑名单挡住
> "白名单扩展名 + 敏感内容"的组合 (如 `kangshan_results.md`), 大小挡住漏网的大二进制。

---

## 2. 白名单 (default-deny: 未列出的一律拒)

**按扩展名 (小写):**

| 类别 | 扩展名 |
|---|---|
| C/C++/CUDA 源码 | `.cu` `.cuh` `.cpp` `.hpp` `.h` `.c` |
| Python | `.py` |
| 构建/配置 | `.cmake` `.toml` `.cfg` `.ini` |
| 脚本 | `.sh` |
| 文档 | `.md`;  `.txt` **仅限仓库根与 `docs/`** (防数据目录里的说明文件夹带) |
| 配置 | `.yml` `.yaml` `.json` (json 仅限**构建/配置**用途 — 程序无法判定语义, 属策略层, 人工把关) |

**按具名文件 (可在任意子目录):** `CMakeLists.txt`、`.gitignore`、`.gitattributes`、`scripts/pre-commit` (钩子自身)。

**按授权路径:** `validation/synthetic/**` —— **合成验证算例**的数据文件放这里。
此目录下文件即便是 `.txt/.csv` 等也按位置放行 (仍受关键词黑名单 + 大小双保险约束)。
**红线: 放进 `validation/synthetic/` 的数据必须是教科书构造 (如解析 dam-break、规则渠道),
无任何真实流域 (鄱阳湖/赣江/…) 来源。** 真实流域数据一律留在仓库外 (`/mnt/f/...`、`~/synxflow_dev/baseline/`)。

---

## 3. 关键词黑名单 (匹配文件路径, 大小写不敏感)

```
poyang 鄱阳  ganjiang 赣江  xinjiang 信江  fuhe 抚河  raohe 饶河
xingzi 星子  duchang 都昌  tangyin  kangshan 康山
case_9   prep_   obs_
```

命中即拦, **即便扩展名在白名单内** (例: `kangshan_note.md` 扩展名合法但命中 `kangshan` → 拦)。

> 取舍说明: 黑名单收录鄱阳湖流域**高辨识度专名** (河流/湖区/测站) 与结构性 token。
> 刻意**不**收 `outlet/thalweg/lean/changjiang/mid/south/east/west` 等通用词, 避免误伤正常源码
> (如 `thalweg` 是通用水力学术语)。这类文件若是真实数据, 其扩展名 (`.asc/.dat/...`) 早已被
> 白名单挡下; default-deny 是兜底。**未发表手稿、核心参数化标定值**等"语义敏感"内容无法靠
> 文件名模式判定, 属人工策略红线 (见 §4), 不要靠钩子兜。

---

## 4. 地形数据红线 & 永久可见性后果

- **地形/测绘数据** (DEM、岸线、断面、`.asc/.nc/.tif/.grd` 栅格) 与**未发表手稿、核心参数化标定**
  绝不进库。这些可能涉及第三方测绘授权/未刊研究, 一旦 GPL+public 即**不可逆**外泄。
- **license 后果**: 本 fork 是 GPL 衍生, push 即公开且不可收回; 真实流域输入/标定数据进库
  = 永久授权他人复制分发。代码可以公开 (这是目标), **数据与未刊成果不行**。
- 钩子只能挡**模式可判定**的部分; 语义敏感内容 (手稿、参数化) 靠本红线 + commit 前自查。

---

## 5. 如何安全放宽白名单

默认拒绝意味着**新增正当文件类型时需要显式改白名单** —— 这是特性不是 bug:
> 误拦十秒可修 (改一行白名单); 漏放 = public + GPL 永久可见的不可逆事故。两者不对称, 故偏向拦。

新增正当类型的步骤:
1. 确认该文件**确无真实流域来源** (源码/构建/通用文档/合成算例);
2. 编辑 `scripts/pre-commit` 的 `ext_ok()` (加扩展名) 或具名/路径分支;
3. 同步更新本文件 §2;
4. 跑一遍 §自测 (见 PR/handoff 里贴的自测结果格式), 确认新类型放行且既有拦截不破;
5. 一并提交 (钩子改动本身受钩子审, `scripts/pre-commit` 在白名单内)。

**不要**为了一次性提交而用 `--no-verify` 绕过 —— 那等于关掉闸门。

---

## 6. 安装 (可复现)

裸 `.git/hooks/` 不进版本库, 故克隆后必须执行一次安装:

```bash
# 推荐 (一行, git 只在该目录找标准 hook 名, 其余脚本忽略):
git config core.hooksPath scripts

# 或软链 (保持 .git/hooks 为标准位置):
ln -sf ../../scripts/pre-commit .git/hooks/pre-commit
chmod +x scripts/pre-commit
```

验证: `git config --get core.hooksPath` 应回 `scripts`; 暂存一个 `test.asc` 跑 `git commit` 应被拦。

---

## 7. 绕过与纪律

- `pre-commit` 是**本地**强约束, `git commit --no-verify` (`-n`) 能绕过它。
- **transport 开发期硬约束: 任何 commit 不得带 `--no-verify` 之类绕过标志。**
  唯一可接受的 "no-" 标志是 `git commit --amend --no-edit` (它不绕过钩子)。
- 服务端无法靠本钩子兜底; 推送前的最后防线就是这道本地闸 + 人工自查 §4。
