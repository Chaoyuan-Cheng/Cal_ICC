# Cal_ICC
# ICC Analysis for Microbiome Data

## Overview

本脚本对物种水平的相对丰度数据进行ICC (Intraclass Correlation Coefficient) 分析，评估Host（个体）、Period（实验阶段）、Diet（饮食）三个因素对肠道菌群组成的影响。

## Input Data

| 文件 | 描述 |
|------|------|
| `relative abundance.csv` | 物种相对丰度矩阵（行=物种，列=样本） |
| `1_metadata-D.txt` | DD组样本元数据 |
| `1_metadata-R.txt` | DR组样本元数据 |

**数据结构**：
- 960个样本（16个个体 x 60天）
- 9140个物种
- 元数据包含：ID, Group, Day, Host, Period, Diet

## Analysis Pipeline

```
1. Load Data
   |
   v
2. Filter Species (mean abundance > 0.5%)
   |  -> 保留24个核心物种
   v
3. For each species:
   |
   |  3.1 Log transform: log10(abundance + 1e-6)
   |
   |  3.2 Fit GAM model:
   |      log_abund ~ s(Day) + s(Host, bs="re") + s(Period, bs="re") + s(Diet, bs="re")
   |
   |  3.3 Extract variance components via gam.vcomp()
   |
   |  3.4 Calculate ICC = Var(effect) / Total_Var
   |
   v
4. Generate Figure
```

## ICC Calculation

ICC衡量某个分组因素能解释多少总变异：

```
ICC = Var(random_effect) / (Var(Host) + Var(Period) + Var(Diet) + Var(residual))
```

- **ICC(Host)**: 个体间差异占总变异的比例
- **ICC(Period)**: 实验阶段差异占总变异的比例
- **ICC(Diet)**: 饮食差异占总变异的比例

## Figure Interpretation

| 元素 | 含义 |
|------|------|
| X轴 | ICC值（0-0.5） |
| Y轴 | 物种名称（按Host ICC降序排列） |
| 三列面板 | Host / Period / Diet 三个效应 |
| 点的颜色 | ICC值大小（黄色=高，蓝色=低） |
| 点的形状 | 显著性（三角形=显著p<0.05，圆形=不显著） |
| 横线长度 | 95%置信区间 |
| 线型 | 实线=显著，虚线=不显著 |
| 灰色填充 | 不显著的效应 |

## Output Files

| 文件 | 描述 |
|------|------|
| `ICC_species.pdf` | 高清矢量图 |
| `ICC_species.png` | PNG格式图片 |
| `ICC_species_results.csv` | ICC数值结果（含置信区间和p值） |

## Parameters

脚本开头可调整的参数：

```r
MIN_MEAN_ABUNDANCE <- 0.005    # 物种筛选阈值 (0.5%)
LINE_WIDTH <- 1.6              # 横线粗细
POINT_SIZE <- 3.0              # 点大小
X_AXIS_LIMIT <- 0.5            # X轴上限
SPECIES_FONT_SIZE <- 7         # 物种名字体大小
```

## Usage

```bash
cd /home/ps/Desktop/work/work_260127/mydata
Rscript --vanilla icc_analysis.R
```
