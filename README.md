# STREAM Python 复现

基于 SoftwareX 论文 10.1016/j.softx.2025.102355 公开模型的 Python 版本。包括 Prototype 和 SoftwareX 的原始 Excel、HiGHS 求解入口、论文情景生成和绘图代码。BSD 版权见 `LICENSE`。

## 首次运行

从 GitHub 下载或克隆后，在 Windows 上安装 64 位 Python 3.12。在本目录打开 PowerShell，执行：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\RUN_PAPER.cmd
```

如果这台电脑已经配置过依赖，也可以使用现有虚拟环境的 Python 来运行本目录的源码。例如：

```powershell
& 'E:\CODE\STREAM-main\STREAM_Python_Reproduction\.venv\Scripts\python.exe' -X utf8 -m paper_reproduction.workflow --only nominal --output new_runs/paper
```

## 计算和输出

- `RUN_PAPER.cmd` 默认从 Excel 重算论文名义情景及其无排放约束的基准情景。每组结果各有 `summary.json`、`annual.csv`、`states.csv`、`trajectory.csv` 和求解日志，位于 `new_runs/paper/runs`。
- `RUN_PAPER.cmd --only all --seconds 300` 求解全部论文情景。`--seconds` 是每个情景的求解时限。结果保留实际状态与最优性差距。
- `RUN_ALL.cmd` 从 Excel 求解原项目的两个案例，结果在 `results`。
- 名义情景算完后，可以执行下列命令生成论文图 3 和图 4。全部情景算完后，省略 `--figures 3 4` 即可绘制图 3 到图 6。

```powershell
.\.venv\Scripts\python.exe -X utf8 -m paper_reproduction.figures --results new_runs/paper --figures 3 4
```

模型和读表代码在 `stream_python`，输入在 `inputs`。`.venv`、计算缓存及历史结果在每台电脑上本地生成，源码包不含这些目录。

## 可选的 Julia 对照和方案改进

`validation/revision_20261110` 含从 Julia 导出的参考矩阵。它们只用于验算，普通 Python 求解时从 Excel 建模。若需要运行 `paper_accuracy_20261111/search_schedules.py`，先将相应情景计算到默认的论文复现历史路径 `paper_reproduction/runs`，例如执行：

```powershell
.\RUN_PAPER.cmd --only nominal --output paper_reproduction
.\.venv\Scripts\python.exe -X utf8 paper_accuracy_20261111/search_schedules.py nominal
```
