# CEL 初始化闭合实验最终报告

## 1. 目的

在已完成轴向边界修正的基准 `F100_G2P20_ZEROPRESSURE_FACEFIX20` 上，只改变 Eulerian 液体初始占据，检验初始流体体积分数是否导致轴向反转。

## 2. 实验身份

- Parent: `F100_G2P20_ZEROPRESSURE_FACEFIX20`
- Child: `F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20`
- 频率 100 Hz，B0 11 mT，G 2.20，S6/S4 ZERO PRESSURE，20 ms
- 网格、几何、接触/摩擦 0.03、材料、阻尼、磁场、输出设置保持不变
- 只执行了一次 Child dynamics run；没有重跑或启动下一候选

## 3. 求解门禁

Child `.sta` 明确写出 `THE ANALYSIS HAS COMPLETED SUCCESSFULLY`，输出 801 个场帧。Datacheck 通过，警告与 Parent 一致：幅值提示、刚性管局部畸变、PRESS 被 S 替代以及接触输出覆盖范围提示。

## 4. 初始液体构造

Parent 使用 7,944 个全液体 Eulerian 单元。Child 使用 8×8×8 中点积分，在名义圆柱腔体内并排除机器人凸包后写入 `*Initial Conditions, type=VOLUME FRACTION`。Child 有 11,517 个正 EVF 单元，其中 3,210 个为部分 EVF。

独立审计得到：实际机器人四面体体积 1.282971394 mm³，凸包体积 1.299769820 mm³；Parent 初始液体 12.186818182 mm³，Child 15.494485973 mm³。按实际四面体几何估计的可达名义液体体积为 15.506303524 mm³，Parent 缺口 21.4073%，Child 缺口 0.0762%。正 EVF 区域均为单一连通分量并轴向连通。

## 5. 几何审计限制

凸包排除对机器人凹陷细节是保守近似；EVF 只规定单元内液体总量，Abaqus 后续界面重构不提供点对点 CAD 布尔保证。16³ 只读复核得到 15.489258229 mm³，与执行的 8³ 结果相差 −0.005227744 mm³。审计中将凸包相对实际四面体造成的 0.016798425 mm³ 过排除单独列出。

## 6. 运动结果

Face-fix Parent 的 C2 在 20 ms 末位移约 −5.5915 mm/s 对应的末端位移为负，反转从约 15.7567 ms 开始，反向位移约 0.03584 mm。Child 的 C1 在约 5.7138 ms 即开始反向；C2 全程为负，20 ms 末速度约 −191.10 mm/s，C2 反向位移约 0.85748 mm，总末端位移约 −0.91366 mm。

近同步 ODB 轨迹比较显示最大位移差约 1.0735 mm、RMS 约 0.4428 mm；最大速度差约 199.3 mm/s、RMS 约 69.4 mm/s。结果属于反转明显放大（结论 C），不是抑制或轻微影响。

## 7. 流体和能量代理量

0–2 ms 湿压绝对值 p95 从约 0.0002845 增至 0.0005041 N/mm²，低端轴向速度代理从约 0.1725 增至 6.6295 mm/s，高端从约 0.0134 增至 0.9912 mm/s。20 ms p95 约 0.01350 N/mm²。Child 20 ms 的 ALLPW 约 0.27625 N·mm、ALLFD 约 0.07625 N·mm、ETOTAL 约 −0.18984 N·mm；这些量只能作为响应证据，不能单独证明失稳或唯一因果。

## 8. 结论

在当前 Face-fix 边界、材料、磁场和接触设置下，修正 Eulerian 初始液体占据没有解决轴向反转，反而显著放大了反转。因此“初始液体占据是唯一主因”被本实验否定。后续若继续隔离，优先做一次接触/罚参数单因素检查；本次交付不启动下一实验。

## 9. 已交付文件

- `initialization_closure/initialization_volume_audit.json`
- `initialization_closure/input_identity_checks.json`
- `initialization_closure/missing_fluid_cell_map.csv`
- `initialization_closure/initialization_closure_summary.json`
- `F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20_robot_motion.gif`
- `F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20_results.zip`

## 10. 辅助日志路径碰撞

Child Fortran 保持原样，仍写入 Parent 目录中的辅助 CSV，因此部分 Parent 辅助日志被 Child 覆盖或追加。Parent ODB 与已导出的 Parent NPZ 未受影响；碰撞文件已归档到 Child 目录并在压缩包清单中标明，不能把这些 CSV 当作未污染的原始 Parent 辅助输出。

## 11. 可复现性

压缩包包含 Child INP、未修改执行 Fortran、磁性表、datacheck/dynamics 文本日志、审计 JSON/CSV、分析脚本和 GIF；约 1 GB 的 ODB、NPZ 原始大文件和 restart 文件不打包，保留在本地结果目录。
