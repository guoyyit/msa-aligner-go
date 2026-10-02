# msa-aligner-go

可复现的多序列比对（MSA）命令行服务。读入一组同源序列（FASTA），输出整体对齐、逐列
一致性分析，并支持并发成对距离计算、确定性合并、崩溃检查点恢复、序列追加与批处理。

## 1. 流程概览

管线由四个确定性阶段组成：

1. **解析（fasta）**：解析并校验 FASTA，规范化大小写与字母表。
2. **距离（distance）**：以可配置并发度计算所有序列对的 k-mer 剖面距离，汇总为
   对称距离矩阵。单对失败（序列过短等）只标记该对，不阻塞其余对。
3. **合并（cluster）**：对距离矩阵做平均连接（UPGMA 风格）凝聚聚类，产生唯一的
   合并顺序。
4. **比对（align）**：严格按合并顺序做带仿射空位的渐进式 profile-profile 比对，
   生成对齐行、一致性序列与逐列统计。

## 2. 输入约束（解析阶段即拒绝）

所有异常输入都会被**明确拒绝**，并返回可区分的错误类型（`fasta.ParseError.Kind`）：

| Kind | 触发条件 |
| --- | --- |
| `missing_id` | `>` 后没有任何标识 |
| `empty_sequence` | 某条记录没有任何合法残基 |
| `illegal_character` | 出现字母表之外的字符（错误信息带行列号与字符） |
| `duplicate_id` | 两条记录首个空白前的标识相同（给出首次出现位置） |
| `malformed` | 首行前出现数据、文件无记录、未知字母表等结构问题 |

字母表：

- `dna`（默认）：`ACGTN-`，允许 `U`（规范为 `T`），小写自动转大写。
- `rna`：`ACGUN-`，允许 `T`（规范为 `U`）。
- `protein`：`ACDEFGHIKLMNPQRSTVWYBXZ*-`。

解析在遇到错误时仍返回已成功解析的记录（部分结果），同时返回聚合错误 `fasta.Errors`，
调用方可以据 `Kind` 区分失败原因。

### 资源上限（在任何重计算前检查）

超限直接返回 `*pipeline.ErrLimitExceeded`，不执行距离/比对：

| 配置项 | 默认值 | 含义 |
| --- | --- | --- |
| `--max-sequences` | `5000` | 最大序列条数 |
| `--max-length` | `200000` | 单条序列最大长度 |
| `--max-total-residues` | `50000000` | 所有残基总数上限 |
| `--max-memory-bytes` | `2147483648` | 距离矩阵 N²×8 字节上限 |

时间与内存可预测：距离为 O(N²) 对、每对 O(L)；矩阵存储 O(N²)；每次 profile 合并为
O(宽₁×宽₂)，宽度受总残基数约束。

## 3. 距离与并发

- 距离定义：`d = 1 - 共同 k-mer 数（取两剖面交集计数）/ min(两序列可用 k-mer 总数)`，
  范围 `[0,1]`，`0` 表示相同。
- DNA/RNA 下 k-mer 按**反向互补规范化**（取 key 与反向互补 key 的较小者），因此
  方向相反的测序读段不会被误判为远距离。
- 不可用对（序列短于 k、无可用 k-mer）矩阵单元记为 `NaN` 并在 `Matrix.Errors` 记录
  原因；这些对不影响其他对计算。
- 并发度由 `--concurrency` 控制（库默认 = CPU 数，CLI 默认 `1` 以保证可观察的顺序）。
  结果只依赖序列内容，**与并发度无关**：任意并发度产生完全相同的矩阵。

## 4. 唯一且可复现的合并顺序

聚类为平均连接：两簇间距离 = 所有有限叶-叶距离按簇大小加权平均。平局与异常处理规则
固定，确保相同输入在任何机器上产生相同合并树：

- 候选按“距离升序，距离相等按 `(较小节点号, 较大节点号)` 升序”选择。
- 没有有限距离的对（NaN）**强制排到最后**；只剩 NaN 候选时，按节点号升序合并，
  记录的距离为 NaN。若全部对均失败则整体拒绝（`cluster: all pairwise distances failed`）。
- 叶节点编号 `0..N-1`，合并产生的新节点按步骤递增编号。

距离矩阵 ⇒ 合并顺序是**纯函数**，因此相同输入必然得到相同合并顺序与相同最终对齐。

## 5. 可配置打分与比对

- 双序列与 profile-profile 均使用仿射空位 Gotoh 三类 DP（对角/上/左）。
- 空位段长度 L 代价 = `gap-open + L * gap-extend`。
- 回溯平局顺序固定为 **对角 > 上 > 左**，保证逐位可复现。
- 默认打分：match `+2`、mismatch `-1`、gap-open `-10`、gap-extend `-1`。
- 参数可配置（`--match/--mismatch/--gap-open/--gap-extend`）；每组参数各自可复现，
  参数改变允许结果改变，但语义一致（同一打分模型）。

## 6. 检查点与崩溃恢复

设置 `--checkpoint <文件>` 后，管线周期性（`--checkpoint-interval`，默认 5s）原子写
检查点：先写临时文件、`fsync`、再 `rename` 覆盖。

文件格式（带完整性尾）：

```
MSACHECKPOINT-V1
<base64(JSON 状态)>
---
<上一行载荷的 hex SHA-256>
```

- 状态内含**输入哈希**（按输入顺序的 ID+数据 SHA-256）与**配置哈希**（字母表、k、
  打分表全枚举、空位罚分的 SHA-256）。
- 重启后若输入或配置不一致，拒绝复用并返回 `ErrCheckpointMismatch`。
- 检查点被截断、魔数错误、base64 损坏、JSON 损坏、**校验和不符**、版本不符，一律
  返回 `checkpoint.ErrCorrupt`，绝不静默恢复。
- 距离阶段按**对**持久化（未完成单元序列化为 JSON `null`，即 NaN）；恢复时只补算
  缺失的对，已完成的对与后续阶段不重复计算。取消（context 取消）也会先落盘已完成部分。
- 聚类比对、最终对齐结果也会缓存；对齐的派生数据（一致性、逐列统计）在恢复时从对齐
  行重新计算，不信任检查点里的派生字段。

## 7. 序列追加（有界且可解释）

`Service.AppendTo` / `align.Append` 在已有对齐上追加新序列：

- **已有行逐字节保持**：从不向已有行插入空位列，残基相对顺序与列位置完全不变。
- 新序列对齐到“冻结的现有一致性 profile”，列**只会新增不会删除或重排**。
- `Alignment.AppendedColumns` 精确记录哪些列是新增的；追加前后宽度差恰好等于该列表
  长度——这就是“与整体重算差异有界且可解释”的保证：差异仅来自新增列，且新增列位置
  被显式给出，不丢失任何已有序列的对齐信息。
- 追加的新 ID 与已有 ID 重复、或新 ID 彼此重复，都会在计算前拒绝；新输入同样经过
  字母表与资源上限校验。

## 8. 进度反馈

进度通过回调（库）或 stderr（CLI）上报，阶段名固定为 `parsing / distance / cluster /
align / done`，并带 `done/total`。批处理还会附带 `batch i/n job=<name>`。

## 9. 输出格式

`--format json`（默认）：

```json
{
  "ids": ["brca1_mouse", "brca1_human"],
  "width": 16,
  "rows": ["ACGTTCGTACGTACGT", "ACGTACGTACGTACGT"],
  "consensus": "ACGT...",
  "columns": [
    {"index": 0, "consensus": "A", "identity": 1.0,
     "conserved": true, "gap_fraction": 0.0, "counts": {"A": 2}}
  ],
  "merge_order": [{"step": 0, "left": 0, "right": 1, "distance": 0.0}],
  "appended_columns": []
}
```

- `columns[].identity`：最常见非空位残基占全部行的比例；并列时取字节序最小的残基。
- `columns[].conserved`：`identity == 1` 且该列无空位。
- `columns[].gap_fraction`：该列空位比例。
- `merge_order`：按执行顺序的合并步骤（`left/right` 为节点号，`distance` 为该步距离，
  NaN 表示该合并不含任何有限对）。

`--format fasta`：仅输出对齐后的 FASTA（`>` 行 + 对齐行）。

## 10. 命令行用法

```bash
# 单个作业（stdin/stdout 用 -）
msa-aligner -input seqs.fasta -concurrency 4 -checkpoint cp.bin -format json -output out.json

# 批处理：目录下所有 .fasta/.fa/.fna，按文件名排序独立执行，
# 输出 <name>.aligned.json（或 .aligned.fasta）
msa-aligner -batch ./jobs -format fasta -concurrency 8

# 崩溃后用相同命令 + 相同 --checkpoint 自动续跑
```

### 全部配置项

| 标志 | 默认值 | 说明 |
| --- | --- | --- |
| `-input` | `-` | 输入 FASTA，`-` 为 stdin |
| `-batch` | 空 | 指定后按批处理目录运行 |
| `-alphabet` | `dna` | `dna` / `rna` / `protein` |
| `-k` | `3` | k-mer 长度 |
| `-concurrency` | `1` | 距离并发 worker 数 |
| `-match` | `2` | 匹配得分 |
| `-mismatch` | `-1` | 错配得分 |
| `-gap-open` | `-10` | 空位开启罚分 |
| `-gap-extend` | `-1` | 空位延伸罚分（每列） |
| `-max-sequences` | `5000` | 序列条数上限 |
| `-max-length` | `200000` | 单序列长度上限 |
| `-max-total-residues` | `50000000` | 总残基上限 |
| `-max-memory-bytes` | `2147483648` | 矩阵内存上限 |
| `-checkpoint` | 空 | 检查点路径（设置即启用恢复） |
| `-checkpoint-interval` | `5s` | 检查点写周期（Go duration） |
| `-output` | `-` | 输出目标，`-` 为 stdout |
| `-format` | `json` | `json` / `fasta` |
| `-quiet` | false | 关闭 stderr 进度 |

退出码：`0` 成功；`1` 作业失败（输入非法、超限、检查点损坏/不匹配等）；`2` 参数/IO 错误。

## 11. 作为库使用

```go
svc := pipeline.New(pipeline.DefaultConfig())
res, err := svc.Align(ctx, f, func(p pipeline.Progress) { log.Printf("%+v", p) })
res.Alignment.Rows      // 对齐行（顺序按合并树成员）
res.Alignment.Columns   // 逐列一致性
res.Tree.Steps          // 唯一合并顺序
```

## 12. 测试

```bash
go test ./...          # 全部测试
go test -race ./...    # 含并发竞态检测
```

覆盖：FASTA 五类异常输入、并发确定性与单对失败隔离、合并平局/NaN 规则、仿射空位
Needleman-Wunsch、渐进式比对、逐列统计、追加不改动已有行、检查点往返与各类损坏检测、
取消后恢复不重复计算、输入/配置不匹配拒绝、参数变化可复现、资源上限、批处理失败隔离。
