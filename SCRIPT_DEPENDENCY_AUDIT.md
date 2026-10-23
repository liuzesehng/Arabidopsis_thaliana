# 脚本自定义检查模块依赖检查

检查范围：当前目录及所有子目录（包括 `trash/`），共 196 个脚本：119 个 Python、65 个 Shell、12 个 R。

已将以下 20 份脚本的外部图形检查改为可选。保留原检查调用和现有检查模块；不新增必须与脚本一起复制的公共模块。

- 缺少 `audit_panel_alignment.py`、`audit_pdf_text.py`、`audit_figure_collisions.py` 或 `validate_figure.py` 时，仅跳过对应检查。
- 终端提示 `[QA SKIPPED]`；对应报告记录 `SKIPPED` 和原因。跳过不会标记为检查通过，汇总报告会标记检查不完整。
- 工具存在时仍执行原检查；执行后的检查失败仍会报错。
- 保留统计计算、数据校验及绘图参数。原始输入数据、常规 Python/R 包及生物信息学工具仍需按原脚本要求准备。

## 修改清单

- [1.transcriptome_analysis/1.10.box.py](1.transcriptome_analysis/1.10.box.py)
- [3.methylome_analysis/3.14.10.bar.contibution.py](3.methylome_analysis/3.14.10.bar.contibution.py)
- [3.methylome_analysis/3.14.11.scatter.box.py](3.methylome_analysis/3.14.11.scatter.box.py)
- [3.methylome_analysis/3.17.tem_meth_plot.py](3.methylome_analysis/3.17.tem_meth_plot.py)
- [3.methylome_analysis/3.19.upset.py](3.methylome_analysis/3.19.upset.py)
- [3.methylome_analysis/3.20.me_tf_plot.py](3.methylome_analysis/3.20.me_tf_plot.py)
- [4.coex_analysis/4.15.TF_CHH_bar.py](4.coex_analysis/4.15.TF_CHH_bar.py)
- [5.trans_experiment/5.6.RCA.scatter.py](5.trans_experiment/5.6.RCA.scatter.py)
- [5.trans_experiment/5.7.RCA.boxplot.py](5.trans_experiment/5.7.RCA.boxplot.py)
- [6.meth_experiment/6.12.scatter.box.py](6.meth_experiment/6.12.scatter.box.py)
- [6.meth_experiment/6.13.tem_meth_plot.py](6.meth_experiment/6.13.tem_meth_plot.py)
- [6.meth_experiment/6.14.meth_tem.py](6.meth_experiment/6.14.meth_tem.py)
- [6.meth_experiment/6.15.meth_tf_loose.py](6.meth_experiment/6.15.meth_tf_loose.py)
- [6.meth_experiment/6.16.chh_tem.py](6.meth_experiment/6.16.chh_tem.py)
- [6.meth_experiment/6.17.chh_filter.py](6.meth_experiment/6.17.chh_filter.py)
- [6.meth_experiment/6.20.chh_tf_plot.py](6.meth_experiment/6.20.chh_tf_plot.py)
- [7.codex_experiment/7.10.chh_trans.py](7.codex_experiment/7.10.chh_trans.py)
- [7.codex_experiment/7.6.module.py](7.codex_experiment/7.6.module.py)
- [7.codex_experiment/7.7.module.contibution.py](7.codex_experiment/7.7.module.contibution.py)
- [7.codex_experiment/7.8.trans_cor.py](7.codex_experiment/7.8.trans_cor.py)

## 完整运行验证（2026-10-10）

按要求清空旧 output 后，已在全新隔离目录重新完成以下验证；本轮未修改脚本源文件，未复用上轮图形或数据表。

20 份涉及自定义图形检查模块的脚本已逐一使用真实输入完整运行，全部通过，共生成 982 个文件（含 78 份 PDF）。验证隐藏外部 nature-figure 检查工具，禁止读取其他自定义 Python 源文件，并将写入限制在对应的隔离目录内。输入文件运行前后 SHA-256 一致；原默认输出路径表达式保持不变。

此前修复：

- 恢复 3.14.11、3.17、3.19、5.6、7.6、7.7 中缺失的 `configure_style` 定义。
- 3.19 跳过对齐检查时，不再搬移未生成的 `alignment.svg`。
- 7.6、7.7 接受缺少检查工具时返回的 `SKIPPED`，不再误判为失败。

119 个 Python、65 个 Shell、12 个 R 脚本的语法检查全部通过；R 使用已有的 `R` conda 环境。其余 176 份脚本完成依赖和语法检查，未运行下载、比对、训练等主体流程。

验证结果位于用户指定的 `list/output/standalone_validation_20261010_195756/`：

- `validation_summary.tsv`：逐份运行结果和输出数量。
- `all_scripts_dependency_scan.tsv`：196 份脚本的依赖与语法检查。
- `output_manifest.tsv`：逐输出文件 SHA-256。
- `original_output_paths_verified.json`：默认输出路径保持不变的证据。
- 各脚本目录中的 `outputs/`、脚本副本和运行日志。

所有生成的 PDF、SVG、PNG、TIFF 均通过读取验证。标准科学计算包和原始输入数据仍需准备；6.14 保留原有 cgmaptools 工具依赖，并已实际执行成功。检查工具缺失记录为 `SKIPPED`，不计为图形 QA 通过。原数据和正式分析结果未覆盖。
