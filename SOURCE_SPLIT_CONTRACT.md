# 生成来源与拆分从属标记

目标：不改变化学输入、计量系数、Flux、节点位置、半径或反应路径，只修正来源与从属的标记。

## 配置

在 `layout.json` 的 `pathway_multiplicity` 中设置：

```json
{"enabled": true, "normalization_flux": 1, "label_mode": "local_flow", "relationship_mode": "source_split"}
```

圆形布局沿用 `instance_marker_style=pool_code`。原来的 `pool` 模式仍可用于复现历史图。

## 三种身份，不可混用

化学身份：代谢物 ID，例如 G3P。重复标签不是新的物理代谢物池。

生成来源：每条实际生成反应的正向计量通量 `S[m,r]*flux[r]/normalization_flux > 0`。
一个来源用 `[A]`；同一代谢物有多个生成反应时，分别用 `[B1]`、`[B2]` 等独立根编号。
编号确定性地按反应 ID 排序，不根据 Candidate 编号或代谢物名字硬编码。

拆分从属：只在唯一已知来源足以覆盖总消耗时，为不含生成的局部使用节点标 `[A.1]`、`[A.2]`。
小数点是父子关系，不是浓度、小数数量或新物质。
主路径直接连接父节点的消耗不增加一个假节点，但必须与远端份额一起计入 consumption_allocations。

因此：`3 F6P [A]` 生成总量中的 1 份沿原主路径使用，其余两个位置各为 `1 F6P [A.1]`、`1 F6P [A.2]`。
`1 G3P [B1]` 和 `1 G3P [B2]` 表示两个不同反应分别生成的贡献，不是从一份 G3P 拆成两份。

## 不确定性和边界

多来源时，单靠稳态计量通量不能确定每个来源精确分配给哪个消耗反应。
消费账本保留 possible_source_ids，但 parent_source_ids 为空，不宣称原子追踪结论。
纯消费的多来源节点使用 `[B:pool]`；实际显示局部消耗大于局部生成时，也使用 `:pool`，
以免把 `max(P,C)` 的消费数错说成某个来源的生成数。
两条生成反应在同一显示节点汇合时可为 `[B1+B2]`；加号表达来源汇合，不表达拆分。
净进口导致来源不足时不强认父子关系；零通量为 `[A:0]`，不虚构来源。
这些标签不表示时间先后、独立浓度池或实测分子轨迹。

## 验证与复现

Python 在原有逐关联分配和守恒校验后推导标记。浏览器从实际 SVG 端点和原始计量/Flux 独立重算；
标签和缓存一起写错仍须失败。数量、来源、从属、来源账本、保存重载、圆形和连通性均接受检查。

```text
python -X utf8 -m unittest test_source_split test_source_split_browser -v
python -X utf8 -m unittest discover -v
```

历史九组完整交付仍保留在本机 `nine_paths_source_split_20260922` 中，
但它属于历史生成结果，不再作为 Git 仓库的运行依赖。新的批量任务应走
根目录稳定入口 `render_batch.py`。只有通过当前仓库验收的结果才可称为本版验收完成。
