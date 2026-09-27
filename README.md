# DNA Assembly API

本地 DNA 序列分析服务（Python / FastAPI）。

- `POST /assemble`：给定 2～12 条读段，求包含全部读段的最短候选序列（最短公共超串）；长度相同时取字典序最小者，结果与输入顺序无关。
- `POST /variants`：给定参考序列与已定位的正向读段，基于碱基质量做对数似然累加，检出高置信度的单碱基替换。
- `POST /align`：给定参考序列与 1～20 条读段，对每条读段及其反向互补做半全局定位，允许替换及至多一次单碱基插入/缺失，返回编辑数不超过 `max_edits` 的最优比对。
- `POST /digest`：给定 12～1000 位序列、线性/环状拓扑与 1～5 个酶切模板（支持 `N/R/Y` 简并碱基），在模板及其反向互补的所有命中处切割，返回去重升序的切点与按 `start` 升序的片段。

## 启动

```bash
docker compose up --build
```

服务监听 `localhost:8000`。

## 接口

### `POST /assemble`

请求体（JSON）：

- `reads`：字符串数组，2～12 条；
- 每条长度 1～80，仅含 `A` / `C` / `G` / `T`（大写）；
- 读段按输入方向使用（不做反向互补），允许重复或相互包含。

```bash
curl -s -X POST http://localhost:8000/assemble \
  -H 'Content-Type: application/json' \
  -d '{"reads": ["ACGT", "GTAC", "TACG"]}'
```

响应：

```json
{
  "sequence": "GTACGT",
  "length": 6,
  "positions": [2, 0, 1]
}
```

- `sequence`：包含全部读段的最短序列；同长度时取字典序最小，且不受输入顺序影响；
- `length`：序列长度；
- `positions`：每条读段在 `sequence` 中最早出现的 0 基下标，按**原输入顺序**返回。

### `POST /variants`

在参考序列上对已定位、无插缺的正向读段逐位做单碱基替换检测。

请求体（JSON）：

- `reference`：1～1000 位，仅含 `A` / `C` / `G` / `T`（大写）；
- `reads`：1～100 条读段；每条包含：
  - `start`：读段在参考序列上的 0 基起始位置（非负）；
  - `sequence`：1～80 位，仅含 `A` / `C` / `G` / `T`；
  - `qualities`：与 `sequence` 等长的整数数组，每个质量值在 2～40；
  - 读段必须完全落在参考范围内（`start + len(sequence) <= len(reference)`）。

打分规则：质量值 Q 对应错读概率 `p = 10^(-Q/10)`；对每个被覆盖位点分别累加四种候选碱基的对数似然——观测碱基与候选相同贡献 `ln(1-p)`，否则贡献 `ln(p/3)`。总分最高的候选胜出（并列按 `A < C < G < T`）；仅当该位点有覆盖、最佳候选不同于参考、且与次高分之差至少 `ln(10)`（10 倍似然比）时才报告变异。

```bash
curl -s -X POST http://localhost:8000/variants \
  -H 'Content-Type: application/json' \
  -d '{
    "reference": "ACGTACGT",
    "reads": [
      {"start": 1, "sequence": "TGT", "qualities": [40, 40, 40]},
      {"start": 1, "sequence": "TGT", "qualities": [38, 39, 40]},
      {"start": 1, "sequence": "TGT", "qualities": [39, 40, 38]}
    ]
  }'
```

响应：

```json
{
  "variants": [
    {"position": 1, "ref": "C", "alt": "T", "depth": 3}
  ]
}
```

- `variants`：按 `position` 升序排列；无覆盖的位点直接忽略，证据不足（差距小于 `ln(10)`）的位点也不返回；
- `position`：变异位点的 0 基坐标；
- `ref` / `alt`：参考碱基 / 检出的替换碱基；
- `depth`：该位点的覆盖读段条数。

### `POST /align`

在参考序列上定位读段：对每条读段及其反向互补，枚举参考的非空连续区间并做全长比对（读段每个碱基都参与比对，不允许末端自由空位）。允许任意替换，外加至多一次单碱基插入或缺失；编辑数 = 替换数 + 插缺数（故插缺至多贡献 1）。

请求体（JSON）：

- `reference`：1～1000 位，仅含 `A` / `C` / `G` / `T`（大写）；
- `reads`：1～20 条读段，每条 1～80 位，仅含 `A` / `C` / `G` / `T`（大写）；
- `max_edits`：整数 0～5，仅返回编辑数不超过它的比对。

```bash
curl -s -X POST http://localhost:8000/align \
  -H 'Content-Type: application/json' \
  -d '{"reference": "ACGTACGT", "reads": ["TACG", "AACGT", "ACGTAA", "CCCC"], "max_edits": 1}'
```

响应：

```json
{
  "alignments": [
    {"start": 1, "end": 5, "strand": "-", "edits": 0, "cigar": "4M"},
    {"start": 0, "end": 5, "strand": "-", "edits": 1, "cigar": "5M"},
    {"start": 0, "end": 6, "strand": "+", "edits": 1, "cigar": "6M"},
    null
  ]
}
```

- `alignments` 按输入顺序返回，每条读段一个元素；无任何达标比对时为 `null`；
- `start` / `end`：命中区间在参考序列上的 0 基半开坐标（`end - start` 为区间长度）；
- `strand`：`+` 表示读段原向匹配，`-` 表示其反向互补匹配（坐标始终在参考正向上）；
- `edits`：替换数加插缺数；
- `cigar`：`M` 为对齐碱基（含替换），`I` 为读段相对参考的插入，`D` 为读段相对参考的缺失；相邻同操作压缩为「计数 + 操作」，如 `2M1I2M`。

最优比对的裁决顺序：编辑数升序 → 插缺数升序 → `start` 升序 → `end` 升序 → 正向（`+`）优先 → `cigar` 字典序（ASCII）最小。

插入 / 缺失示例（读段 `AACGC` 在区间 `[0,4)` 插入一个 `C`；读段 `AGACG` 在区间 `[1,7)` 跳过一个参考碱基）：

```bash
curl -s -X POST http://localhost:8000/align \
  -H 'Content-Type: application/json' \
  -d '{"reference": "AAGCACG", "reads": ["AACGC", "AGACG"], "max_edits": 1}'
```

```json
{
  "alignments": [
    {"start": 0, "end": 4, "strand": "+", "edits": 1, "cigar": "2M1I2M"},
    {"start": 1, "end": 7, "strand": "+", "edits": 1, "cigar": "2M1D3M"}
  ]
}
```

### `POST /digest`

限制性内切酶式的电子酶切：每个模板描述一个识别位点及其切点偏移，模板本身与其**反向互补**的每次命中都会产生一个切点。

请求体（JSON）：

- `sequence`：12～1000 位，仅含 `A` / `C` / `G` / `T`（大写）；
- `topology`：`"linear"`（线性）或 `"circular"`（环状）；
- `templates`：1～5 项，每项包含：
  - `pattern`：2～12 位，字符取自 `A` / `C` / `G` / `T` / `N` / `R` / `Y`（大写）；`N` 匹配任意碱基，`R` 匹配 `A/G`，`Y` 匹配 `C/T`；
  - `offset`：0～`len(pattern)` 的整数。

切点规则（命中起点为 0 基坐标 `s`）：

- 模板正向命中：切点为 `s + offset`；
- 模板反向互补命中：切点为 `s + len(pattern) - offset`；
- 线性序列的命中必须完整落在序列内，切点范围为 `0..len(sequence)`（落在端点 0 / 末端的切点保留在 `cuts` 中，但不产生空片段）；
- 环状序列的命中**可以跨过序列末尾**（把序列视为首尾相接），切点对序列长度取模。

线性示例（EcoRI 风格的回文位点 `GAATTC`，`offset=1`）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAGAATTCTTAA",
    "topology": "linear",
    "templates": [{"pattern": "GAATTC", "offset": 1}]
  }'
```

响应：

```json
{
  "cuts": [3, 7],
  "fragments": [
    {"start": 0, "end": 3, "sequence": "AAG"},
    {"start": 3, "end": 7, "sequence": "AATT"},
    {"start": 7, "end": 12, "sequence": "CTTAA"}
  ]
}
```

环状示例（回文位点 `GTAC` 在环上出现两次，`offset=2` 为正中切割；末段跨零点）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAAAGTACGTAC",
    "topology": "circular",
    "templates": [{"pattern": "GTAC", "offset": 2}]
  }'
```

```json
{
  "cuts": [6, 10],
  "fragments": [
    {"start": 6, "end": 10, "sequence": "ACGT"},
    {"start": 10, "end": 6, "sequence": "ACAAAAGT"}
  ]
}
```

环状序列无切点时，返回一个完整片段 `(start=0, end=0, 原序列)`；恰好一个切点时，返回从该切点起读完整圈的一个片段（`start == end == 切点`）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "CCCCCCCCCCCC",
    "topology": "circular",
    "templates": [{"pattern": "AATT", "offset": 1}]
  }'
```

```json
{
  "cuts": [],
  "fragments": [
    {"start": 0, "end": 0, "sequence": "CCCCCCCCCCCC"}
  ]
}
```

响应字段：

- `cuts`：所有模板正向与反向互补命中产生的切点，**去重后按数值升序**；环状为 `0..len(sequence)-1`（取模），线性为 `0..len(sequence)`；
- `fragments`：按 `start` 升序，每项含 `start` / `end`（0 基半开坐标）与 `sequence`：
  - 线性：以 `0`、各切点、序列末端为界分段，落在两端的切点只会产生被丢弃的空片段，故无空片段；
  - 环状：相邻切点之间各成一段，最后一段跨过零点（此时 `end < start`，序列按环绕顺序拼接）；无切点与单切点按上述特例返回。

### 错误处理

- `/assemble`：读段数量不在 2～12、长度不在 1～80，或含 `A/C/G/T` 以外字符时，返回 4xx（422）及错误详情。
- `/variants`：参考序列或读段长度越界、含非法碱基、质量值越界、`sequence` 与 `qualities` 不等长、读段超出参考范围、`start` 为负等，均返回 4xx（422）及错误详情。
- `/align`：参考序列长度不在 1～1000、读段数量不在 1～20、读段长度不在 1～80、含 `A/C/G/T` 以外字符（含小写）、`max_edits` 不是 0～5 的整数，或缺字段、类型错误等，均返回 4xx（422）及错误详情；读段无达标比对不报错，对应结果为 `null`。
- `/digest`：`sequence` 长度不在 12～1000 或含 `A/C/G/T` 以外字符（含小写）、`topology` 不是 `linear`/`circular`、`templates` 数量不在 1～5、`pattern` 长度不在 2～12 或含 `A/C/G/T/N/R/Y` 以外字符、`offset` 不是 0～`len(pattern)` 的整数、`pattern` 长于 `sequence`，或缺字段、类型错误等，均返回 4xx（422）及错误详情；没有任何命中不报错，按无切点规则返回。

### 其他

- `GET /health`：健康检查，返回 `{"status": "ok"}`；
- `GET /docs`：FastAPI 自动生成的交互式接口文档。

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload
```

运行测试：

```bash
pip install pytest httpx
python -m pytest tests/ -v
```

## 实现说明

- 先去除完全重复的读段，再移除被其他读段包含的读段（不影响超串结果）；
- 计算读段两两之间的最大前后缀重叠，用子集 DP（Held-Karp 风格，状态数 ≤ 2^12 × 12）枚举所有拼接顺序，按「长度最短，其次字典序最小」取最优；
- 字典序按 ASCII（`A < C < G < T`），因此结果对输入顺序、重复读段均不变。
- `/variants` 的计算位于独立模块 `app/variants.py`，不依赖 Web 层；`ln(1-p)` 与 `ln(p/3)` 对质量值 2～40 预先查表，再按覆盖关系逐位累加四个候选的对数似然；遍历候选时严格使用「大于」更新最佳下标，天然实现 `A < C < G < T` 的并列裁决。
- `/align` 的计算位于独立模块 `app/aligner.py`，不依赖 Web 层：对正向与反向互补两条序列，分别枚举每个起点上的无空位比对（替换）与一次单碱基插入/缺失比对，插入/缺失两侧的错配数用左右两个错配计数数组线性累加到 `O(读长)`；只保留 `编辑数 ≤ max_edits` 的候选，并按 `(edits, indels, start, end, strand, cigar)` 取最小。剪枝包括：错配计数达到 `max_edits` 后该方向停止扩展、候选在 `(edits, indels, start, end)` 前缀上已落败则不构造 CIGAR、找到起点最小的完美无空位比对后直接结束该链搜索。
- `/digest` 的计算位于独立模块 `app/digester.py`，不依赖 Web 层：把序列复制一份（`sequence + sequence`）即可让环状命中与跨零点片段直接按下标读取，无需逐位取模；简并字符用 `N→ACGT / R→AG / Y→CT` 的字符集合判定，反向互补时 `N→N`、`R↔Y`。每个模板在每个候选起点分别尝试正向模式与其反向互补，切点汇入一个集合天然去重。线性分段以 `0、切点、末端` 为界并丢弃零长度段；环状无切点返回 `(0, 0, 原序列)`，单切点返回该点起的整圈，多切点在相邻切点间分段（末段 `end < start`，跨零点）。
