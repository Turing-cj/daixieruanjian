# 局部数量标注 Contract

本规则修正共享代谢物在不同绘图位置重复显示整池总量的问题。
它只处理标注和可审计元数据，不改模型、反应计量、Flux、节点位置、连线或主圆。

## 启用

在标准 `layout.json` 中配置：

```json
"pathway_multiplicity": {
  "enabled": true,
  "normalization_flux": 1,
  "label_mode": "local_flow"
}
```

之后仍使用原来的 `render_network.py` 入口。
旧布局未指定 `label_mode` 时保留 `pool_total` 行为，便于复现历史结果；新图使用 `local_flow`。
本规则没有 Candidate 编号、F6P、G3P 或特定反应名称的条件分支。

## 数值语义

对代谢物 m、反应 r，归一化的有向计量通量为 `S[m,r] * flux[r] / normalization_flux`。
正值是生成，负值是消耗。反向通量、非单位计量系数、小数和零都按同一规则处理。

每个反应-代谢物关联只归入一个实际显示节点：普通节点、局部副本或共享局部节点。
共享局部节点先合并它连接的反应，再分别求总生成量 P、总消耗量 C；标签为 `max(P,C)`。
不是 `P+C`，否则同一份物质从生成到消耗会被重复计算。
若一个上游位置生成了总量 3，而另一个局部位置只消耗 1，它们分别显示 3 与 1。
相同 `[A]` 等标记仍引用同一个化学代谢物池，不是彼此独立的新池。

`pathway_multiplicity.counts` 和 `balances` 保留全池总量、总生成/消耗和净值。
`view_flows` 保留每个显示节点的局部数量、全池总量及逐反应有向计量通量。
每个显示节点还有 `pool_total` 与 `amount_scope=local_flow` 元数据。
生产侧与消费侧应分别加总回全池对应侧；不可将所有标签数字直接相加。
这些数字是归一化通量计量数，不是浓度、独立分子池或绝对实测通量。

## 验收

Python 校验每个关联只分配一次、各侧加总守恒、共享副本身份一致，全部验证后才写入标注。
浏览器不信任缓存：从实际 SVG 连线端点和原始计量/Flux 独立重算局部数量。
同时保留圆形、箭头、连通性、断线、交叉、标签碰撞和越界验收。
错误标注即使与缓存同步修改也不得通过。

## 当前仓库验收

局部数量语义由根目录回归测试直接覆盖，当前支持的验证入口是：

```text
python -X utf8 -m unittest test_branch_allocation test_branch_browser -v
python -X utf8 -m unittest discover -v
```

历史九组完整交付仍保留在本机
`nine_paths_branch_amounts_20260922` 目录中，但它属于历史生成结果，不再作为
Git 仓库的运行依赖。新的批量任务应走根目录稳定入口 `render_batch.py`。
