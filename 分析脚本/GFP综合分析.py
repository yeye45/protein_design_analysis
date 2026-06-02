from __future__ import annotations

import html
import math
import re
import shutil
from copy import copy
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


根目录 = Path(__file__).resolve().parents[1]
输入文件 = 根目录 / "GFP_data.xlsx"
参考序列文件 = 根目录 / "AAseqs of 5 GFP proteins.txt"
输出目录 = 根目录 / "分析结果"
图表目录 = 输出目录 / "图表"
工作簿路径 = 输出目录 / "GFP综合分析结果.xlsx"
报告路径 = 输出目录 / "GFP综合分析报告.html"

类型顺序 = ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]
类型颜色 = {
    "avGFP": "#2A9D8F",
    "amacGFP": "#457B9D",
    "cgreGFP": "#E76F51",
    "ppluGFP": "#8A5A44",
}
氨基酸顺序 = list("ACDEFGHIKLMNPQRSTVWY")
突变格式 = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
检测下限阈值 = 1.30104


def 设置绘图样式() -> None:
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 180
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False


def 读取参考序列() -> dict[str, str]:
    序列: dict[str, str] = {}
    当前名称: str | None = None
    for 原始行 in 参考序列文件.read_text(encoding="utf-8").splitlines():
        行 = 原始行.strip()
        if not 行:
            continue
        if 行.startswith(">"):
            当前名称 = 行[1:]
            序列[当前名称] = ""
        elif 当前名称:
            序列[当前名称] += 行
    return 序列


def 拆分突变(突变组合: str) -> list[str]:
    return [] if 突变组合 == "WT" else str(突变组合).split(":")


def 保存图片(文件名: str) -> None:
    plt.tight_layout()
    plt.savefig(图表目录 / 文件名, bbox_inches="tight")
    plt.close()


def 百分比(数值: float) -> str:
    return f"{数值:.2%}"


def 列表转文本(值: list[str]) -> str:
    return ":".join(值) if 值 else "WT"


def 构造基础数据() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str], pd.DataFrame]:
    亮度 = pd.read_excel(输入文件, sheet_name="brightness")
    候选 = pd.read_excel(输入文件, sheet_name="beforetopseqs")
    参考序列 = 读取参考序列()

    亮度["突变数量"] = 亮度["aaMutations"].map(lambda x: len(拆分突变(x)))
    wt = (
        亮度.loc[亮度["aaMutations"].eq("WT"), ["GFP type", "Brightness"]]
        .set_index("GFP type")["Brightness"]
        .to_dict()
    )
    亮度["WT亮度"] = 亮度["GFP type"].map(wt)
    亮度["相对WT亮度变化"] = 亮度["Brightness"] - 亮度["WT亮度"]
    亮度["疑似检测下限"] = 亮度["Brightness"].le(检测下限阈值)
    亮度["高于WT"] = 亮度["相对WT亮度变化"].gt(0)

    数量基线 = (
        亮度.groupby(["GFP type", "突变数量"], as_index=False)["Brightness"]
        .mean()
        .rename(columns={"Brightness": "同突变数量平均亮度"})
    )
    亮度 = 亮度.merge(数量基线, on=["GFP type", "突变数量"], how="left")
    亮度["校正后亮度残差"] = 亮度["Brightness"] - 亮度["同突变数量平均亮度"]

    展开记录: list[dict[str, object]] = []
    for 行号, 行 in 亮度.iterrows():
        for token in 拆分突变(行["aaMutations"]):
            匹配 = 突变格式.match(token)
            if not 匹配:
                原始氨基酸, 数据位点, 新氨基酸 = "?", np.nan, "?"
            else:
                原始氨基酸, 位点文本, 新氨基酸 = 匹配.groups()
                数据位点 = int(位点文本)
            展开记录.append(
                {
                    "原记录行号": 行号 + 2,
                    "GFP类型": 行["GFP type"],
                    "突变组合": 行["aaMutations"],
                    "突变": token,
                    "原始氨基酸": 原始氨基酸,
                    "数据位点_0based": 数据位点,
                    "常规位点_1based": 数据位点 + 1 if pd.notna(数据位点) else np.nan,
                    "新氨基酸": 新氨基酸,
                    "突变数量": 行["突变数量"],
                    "亮度": 行["Brightness"],
                    "相对WT亮度变化": 行["相对WT亮度变化"],
                    "校正后亮度残差": 行["校正后亮度残差"],
                    "高于WT": 行["高于WT"],
                }
            )
    展开 = pd.DataFrame(展开记录)
    return 亮度, 候选, 参考序列, 展开


def 生成数据质量表(
    亮度: pd.DataFrame, 候选: pd.DataFrame, 参考序列: dict[str, str], 展开: pd.DataFrame
) -> pd.DataFrame:
    质量记录: list[dict[str, object]] = []

    def 添加(检查项: str, 结果: object, 说明: str) -> None:
        质量记录.append({"检查项": 检查项, "结果": 结果, "说明": 说明})

    添加("brightness 数据行数", len(亮度), "不含表头")
    添加("beforetopseqs 数据行数", len(候选), "不含表头")
    添加("brightness 缺失值数量", int(亮度[["aaMutations", "GFP type", "Brightness"]].isna().sum().sum()), "三个原始字段")
    添加("beforetopseqs 缺失值数量", int(候选[["year", "sequence"]].isna().sum().sum()), "两个原始字段")
    添加("brightness 完全重复行", int(亮度[["aaMutations", "GFP type", "Brightness"]].duplicated().sum()), "预期为 0")
    添加("类型内重复突变组合", int(亮度[["GFP type", "aaMutations"]].duplicated().sum()), "预期为 0")
    添加("候选序列完全重复行", int(候选.duplicated().sum()), "预期为 0")
    添加("疑似检测下限记录", int(亮度["疑似检测下限"].sum()), f"亮度 <= {检测下限阈值}")
    添加("最大突变数量", int(亮度["突变数量"].max()), "用于识别极端多点突变")
    添加("解析失败的突变 token", int((展开["原始氨基酸"] == "?").sum()), "支持普通突变及 *238G 一类终止位点替换")
    添加("终止位点替换 token", int((展开["原始氨基酸"] == "*").sum()), "保留为特殊突变")

    普通突变 = 展开.loc[展开["原始氨基酸"].ne("*") & 展开["原始氨基酸"].ne("?")].copy()
    普通突变["参考序列原始氨基酸"] = 普通突变.apply(
        lambda x: 参考序列[x["GFP类型"]][int(x["数据位点_0based"])], axis=1
    )
    添加(
        "按 0-based 编号无法匹配参考序列的普通 token",
        int((普通突变["原始氨基酸"] != 普通突变["参考序列原始氨基酸"]).sum()),
        "数据位点 0 对应参考序列第 1 位；预期为 0",
    )
    添加("候选序列长度集合", ", ".join(map(str, sorted(候选["sequence"].str.len().unique()))), "候选序列应保持相同长度")
    非法字符 = sorted(set("".join(候选["sequence"])) - set(氨基酸顺序))
    添加("候选序列非法字符", ", ".join(非法字符) if 非法字符 else "无", "仅允许标准 20 种氨基酸")
    for 类型 in 类型顺序:
        添加(f"{类型} 参考序列长度", len(参考序列[类型]), "来自 AAseqs of 5 GFP proteins.txt")
    return pd.DataFrame(质量记录)


def 生成类型概览(亮度: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    概览行 = []
    for 类型 in 类型顺序:
        数据 = 亮度.loc[亮度["GFP type"].eq(类型)]
        wt = 数据.loc[数据["aaMutations"].eq("WT"), "Brightness"].iloc[0]
        概览行.append(
            {
                "GFP类型": 类型,
                "样本数量": len(数据),
                "WT亮度": wt,
                "最小亮度": 数据["Brightness"].min(),
                "最大亮度": 数据["Brightness"].max(),
                "亮度范围": 数据["Brightness"].max() - 数据["Brightness"].min(),
                "平均亮度": 数据["Brightness"].mean(),
                "亮度中位数": 数据["Brightness"].median(),
                "亮度标准差": 数据["Brightness"].std(),
                "疑似检测下限数量": int(数据["疑似检测下限"].sum()),
                "疑似检测下限比例": 数据["疑似检测下限"].mean(),
                "高于WT数量": int(数据["高于WT"].sum()),
                "高于WT比例": 数据["高于WT"].mean(),
                "突变数量中位数": 数据["突变数量"].median(),
                "突变数量最大值": 数据["突变数量"].max(),
            }
        )
    概览 = pd.DataFrame(概览行)
    分位数 = (
        亮度.groupby("GFP type")["Brightness"]
        .quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1])
        .unstack()
        .reset_index()
        .rename(
            columns={
                "GFP type": "GFP类型",
                0: "最小值",
                0.01: "P01",
                0.05: "P05",
                0.25: "P25",
                0.5: "P50",
                0.75: "P75",
                0.95: "P95",
                0.99: "P99",
                1: "最大值",
            }
        )
    )
    数量统计 = (
        亮度.groupby(["GFP type", "突变数量"])["Brightness"]
        .agg(["size", "mean", "median", "std", "min", "max"])
        .reset_index()
        .rename(
            columns={
                "GFP type": "GFP类型",
                "size": "样本数量",
                "mean": "平均亮度",
                "median": "亮度中位数",
                "std": "亮度标准差",
                "min": "最小亮度",
                "max": "最大亮度",
            }
        )
    )
    return 概览, 分位数, 数量统计


def 生成替换与位点关联表(展开: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    替换 = (
        展开.groupby(["GFP类型", "突变", "原始氨基酸", "数据位点_0based", "常规位点_1based", "新氨基酸"])
        .agg(
            出现次数=("突变", "size"),
            平均亮度=("亮度", "mean"),
            平均相对WT变化=("相对WT亮度变化", "mean"),
            平均校正后残差=("校正后亮度残差", "mean"),
            校正后残差中位数=("校正后亮度残差", "median"),
            校正后残差标准差=("校正后亮度残差", "std"),
            高于WT比例=("高于WT", "mean"),
        )
        .reset_index()
    )
    替换["价值可靠性"] = np.where(替换["出现次数"].ge(10), "较可靠", "样本偏少")
    替换["价值分类"] = np.select(
        [
            替换["出现次数"].ge(10) & 替换["平均校正后残差"].ge(0.05),
            替换["出现次数"].ge(10) & 替换["平均校正后残差"].le(-0.05),
        ],
        ["高价值", "低价值"],
        default="中性或证据不足",
    )
    替换 = 替换.sort_values(["GFP类型", "平均校正后残差"], ascending=[True, False])

    位点 = (
        展开.groupby(["GFP类型", "数据位点_0based", "常规位点_1based"])
        .agg(
            出现次数=("突变", "size"),
            涉及替换种类=("突变", "nunique"),
            平均亮度=("亮度", "mean"),
            平均相对WT变化=("相对WT亮度变化", "mean"),
            平均校正后残差=("校正后亮度残差", "mean"),
            校正后残差中位数=("校正后亮度残差", "median"),
            高于WT比例=("高于WT", "mean"),
        )
        .reset_index()
    )
    位点["价值可靠性"] = np.where(位点["出现次数"].ge(30), "较可靠", "样本偏少")
    位点["价值分类"] = np.select(
        [
            位点["出现次数"].ge(30) & 位点["平均校正后残差"].ge(0.03),
            位点["出现次数"].ge(30) & 位点["平均校正后残差"].le(-0.03),
        ],
        ["高价值", "低价值"],
        default="中性或证据不足",
    )
    位点 = 位点.sort_values(["GFP类型", "平均校正后残差"], ascending=[True, False])

    替换类型 = (
        展开.assign(替换类型=展开["原始氨基酸"] + ">" + 展开["新氨基酸"])
        .groupby(["GFP类型", "替换类型"])
        .agg(
            出现次数=("突变", "size"),
            平均相对WT变化=("相对WT亮度变化", "mean"),
            平均校正后残差=("校正后亮度残差", "mean"),
            高于WT比例=("高于WT", "mean"),
        )
        .reset_index()
        .sort_values(["GFP类型", "平均校正后残差"], ascending=[True, False])
    )
    return 替换, 位点, 替换类型


def 生成全类型价值表(
    位点关联: pd.DataFrame,
    替换类型: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    可靠位点 = 位点关联.loc[位点关联["出现次数"].ge(30)].copy()
    位点透视 = (
        可靠位点.pivot_table(
            index="常规位点_1based",
            columns="GFP类型",
            values="平均校正后残差",
            aggfunc="mean",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    位点汇总 = (
        可靠位点.groupby("常规位点_1based")
        .agg(
            覆盖GFP类型数=("GFP类型", "nunique"),
            覆盖GFP类型=("GFP类型", lambda x: ",".join(sorted(set(x)))),
            总出现次数=("出现次数", "sum"),
            平均校正后残差_类型等权=("平均校正后残差", "mean"),
            中位校正后残差_类型等权=("平均校正后残差", "median"),
            高价值类型数=("平均校正后残差", lambda x: int((x >= 0.03).sum())),
            低价值类型数=("平均校正后残差", lambda x: int((x <= -0.03).sum())),
        )
        .reset_index()
    )
    位点汇总 = 位点汇总.merge(位点透视, on="常规位点_1based", how="left")
    位点汇总["跨类型价值分类"] = np.select(
        [
            位点汇总["覆盖GFP类型数"].ge(2)
            & 位点汇总["平均校正后残差_类型等权"].ge(0.03)
            & 位点汇总["高价值类型数"].ge(2)
            & 位点汇总["低价值类型数"].eq(0),
            位点汇总["覆盖GFP类型数"].ge(2)
            & 位点汇总["平均校正后残差_类型等权"].le(-0.03)
            & 位点汇总["低价值类型数"].ge(2)
            & 位点汇总["高价值类型数"].eq(0),
            位点汇总["覆盖GFP类型数"].ge(2)
            & 位点汇总["高价值类型数"].gt(0)
            & 位点汇总["低价值类型数"].gt(0),
        ],
        ["跨类型高价值", "跨类型低价值", "跨类型分歧"],
        default="证据不足或中性",
    )
    位点汇总 = 位点汇总.sort_values("平均校正后残差_类型等权", ascending=False)

    替换类型汇总 = (
        替换类型.groupby("替换类型")
        .agg(
            覆盖GFP类型数=("GFP类型", "nunique"),
            覆盖GFP类型=("GFP类型", lambda x: ",".join(sorted(set(x)))),
            总出现次数=("出现次数", "sum"),
            平均相对WT变化_类型等权=("平均相对WT变化", "mean"),
            平均校正后残差_类型等权=("平均校正后残差", "mean"),
            中位校正后残差_类型等权=("平均校正后残差", "median"),
            高价值类型数=("平均校正后残差", lambda x: int((x >= 0.03).sum())),
            低价值类型数=("平均校正后残差", lambda x: int((x <= -0.03).sum())),
        )
        .reset_index()
    )
    替换类型汇总["跨类型价值分类"] = np.select(
        [
            替换类型汇总["覆盖GFP类型数"].ge(2)
            & 替换类型汇总["总出现次数"].ge(20)
            & 替换类型汇总["平均校正后残差_类型等权"].ge(0.03)
            & 替换类型汇总["高价值类型数"].ge(2)
            & 替换类型汇总["低价值类型数"].eq(0),
            替换类型汇总["覆盖GFP类型数"].ge(2)
            & 替换类型汇总["总出现次数"].ge(20)
            & 替换类型汇总["平均校正后残差_类型等权"].le(-0.03)
            & 替换类型汇总["低价值类型数"].ge(2)
            & 替换类型汇总["高价值类型数"].eq(0),
            替换类型汇总["覆盖GFP类型数"].ge(2)
            & 替换类型汇总["高价值类型数"].gt(0)
            & 替换类型汇总["低价值类型数"].gt(0),
        ],
        ["跨类型高价值", "跨类型低价值", "跨类型分歧"],
        default="证据不足或中性",
    )
    替换类型汇总 = 替换类型汇总.sort_values("平均校正后残差_类型等权", ascending=False)

    return 位点汇总, 替换类型汇总


def 生成共现组合表(亮度: pd.DataFrame) -> pd.DataFrame:
    汇总: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for _, 行 in 亮度.iterrows():
        tokens = 拆分突变(行["aaMutations"])
        for 左, 右 in combinations(sorted(tokens), 2):
            汇总[(行["GFP type"], 左, 右)].append(float(行["校正后亮度残差"]))
    记录 = []
    for (类型, 左, 右), 残差列表 in 汇总.items():
        if len(残差列表) < 8:
            continue
        记录.append(
            {
                "GFP类型": 类型,
                "突变1": 左,
                "突变2": 右,
                "共同出现次数": len(残差列表),
                "组合所在序列平均校正后残差": np.mean(残差列表),
                "组合所在序列残差中位数": np.median(残差列表),
            }
        )
    if not 记录:
        return pd.DataFrame(columns=["GFP类型", "突变1", "突变2", "共同出现次数", "组合所在序列平均校正后残差", "组合所在序列残差中位数"])
    return pd.DataFrame(记录).sort_values(
        ["GFP类型", "组合所在序列平均校正后残差"], ascending=[True, False]
    )


def 序列差异(参考: str, 序列: str) -> list[tuple[int, str, str]]:
    return [(i + 1, 原始, 新) for i, (原始, 新) in enumerate(zip(参考, 序列)) if 原始 != 新]


def 生成候选序列分析(
    候选: pd.DataFrame,
    参考序列: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    av参考 = 参考序列["avGFP"]
    临时参考 = 候选.iloc[0]["sequence"]
    候选记录 = []
    明细记录 = []
    for 序号, 行 in 候选.reset_index(drop=True).iterrows():
        序列 = 行["sequence"]
        av差异 = 序列差异(av参考, 序列)
        临时参考差异 = 序列差异(临时参考, 序列)
        tokens = [f"{原始}{位点 - 1}{新}" for 位点, 原始, 新 in av差异]
        候选记录.append(
            {
                "候选序号": 序号 + 1,
                "年份": 行["year"],
                "序列长度": len(序列),
                "相对avGFP参考突变数": len(av差异),
                "相对avGFP参考突变": 列表转文本(tokens),
                "相对2024首条序列差异数": len(临时参考差异),
                "相对2024首条序列差异": 列表转文本([f"{原始}{位点}{新}" for 位点, 原始, 新 in 临时参考差异]),
                "序列": 序列,
            }
        )
        for 位点, 原始, 新 in av差异:
            token = f"{原始}{位点 - 1}{新}"
            明细记录.append(
                {
                    "候选序号": 序号 + 1,
                    "年份": 行["year"],
                    "数据突变编号_0based": token,
                    "常规位点_1based": 位点,
                    "原始氨基酸": 原始,
                    "新氨基酸": 新,
                }
            )
    候选表 = pd.DataFrame(候选记录).sort_values(["年份", "候选序号"])
    明细表 = pd.DataFrame(明细记录)
    if 明细表.empty:
        年度频率 = pd.DataFrame(columns=["年份", "数据突变编号_0based", "出现序列数", "年度序列数", "年度频率"])
        年度变化 = pd.DataFrame(columns=["数据突变编号_0based", "2024频率", "2025频率", "频率变化_2025减2024"])
    else:
        年度总数 = 候选.groupby("year").size().to_dict()
        年度频率 = (
            明细表.groupby(["年份", "数据突变编号_0based"])
            .size()
            .reset_index(name="出现序列数")
        )
        年度频率["年度序列数"] = 年度频率["年份"].map(年度总数)
        年度频率["年度频率"] = 年度频率["出现序列数"] / 年度频率["年度序列数"]
        年度变化 = (
            年度频率.pivot(index="数据突变编号_0based", columns="年份", values="年度频率")
            .fillna(0)
            .reset_index()
            .rename(columns={2024: "2024频率", 2025: "2025频率"})
        )
        for 列 in ["2024频率", "2025频率"]:
            if 列 not in 年度变化:
                年度变化[列] = 0.0
        年度变化["频率变化_2025减2024"] = 年度变化["2025频率"] - 年度变化["2024频率"]
        年度变化 = 年度变化.sort_values("频率变化_2025减2024", ascending=False)
    return 候选表, 明细表, 年度频率, 年度变化


def 绘制类型亮度范围(亮度: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    数据列表 = [亮度.loc[亮度["GFP type"].eq(x), "Brightness"] for x in 类型顺序]
    bp = ax.boxplot(数据列表, vert=False, tick_labels=类型顺序, patch_artist=True, showfliers=False)
    for 箱, 类型 in zip(bp["boxes"], 类型顺序):
        箱.set_facecolor(类型颜色[类型])
        箱.set_alpha(0.65)
    for y, 类型 in enumerate(类型顺序, start=1):
        数据 = 亮度.loc[亮度["GFP type"].eq(类型)]
        wt = 数据.loc[数据["aaMutations"].eq("WT"), "Brightness"].iloc[0]
        ax.scatter([wt], [y], marker="*", s=180, color="#222222", zorder=4)
        ax.plot([数据["Brightness"].min(), 数据["Brightness"].max()], [y, y], color="#555555", lw=1, alpha=0.7)
    ax.set_title("各 GFP 类型亮度范围与分布")
    ax.set_xlabel("亮度")
    ax.set_ylabel("GFP 类型")
    ax.grid(axis="x", alpha=0.25)
    ax.text(0.995, 0.03, "黑色星标：WT 亮度", transform=ax.transAxes, ha="right", fontsize=9)
    保存图片("各类型亮度范围图.png")


def 绘制类型亮度分布(亮度: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False)
    for ax, 类型 in zip(axes.flat, 类型顺序):
        数据 = 亮度.loc[亮度["GFP type"].eq(类型), "Brightness"]
        wt = 亮度.loc[亮度["GFP type"].eq(类型) & 亮度["aaMutations"].eq("WT"), "Brightness"].iloc[0]
        ax.hist(数据, bins=60, color=类型颜色[类型], alpha=0.78, edgecolor="white", linewidth=0.3)
        ax.axvline(wt, color="#222222", linestyle="--", linewidth=1.5, label=f"WT={wt:.3f}")
        ax.set_title(f"{类型}（n={len(数据):,}）")
        ax.set_xlabel("亮度")
        ax.set_ylabel("样本数")
        ax.legend(frameon=False)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("各 GFP 类型亮度直方图", fontsize=14)
    保存图片("各类型亮度分布图.png")


def 绘制突变数量关系(数量统计: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for 类型 in 类型顺序:
        数据 = 数量统计.loc[数量统计["GFP类型"].eq(类型)].sort_values("突变数量")
        axes[0].plot(数据["突变数量"], 数据["平均亮度"], marker="o", ms=3.5, label=类型, color=类型颜色[类型])
        axes[1].plot(数据["突变数量"], 数据["样本数量"], marker="o", ms=3.5, label=类型, color=类型颜色[类型])
    axes[0].set_ylabel("平均亮度")
    axes[0].set_title("突变数量与平均亮度关系")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, ncol=4)
    axes[1].set_xlabel("突变数量")
    axes[1].set_ylabel("样本数（对数轴）")
    axes[1].set_yscale("log")
    axes[1].grid(alpha=0.25)
    保存图片("突变数量与亮度关系图.png")


def 绘制位点价值(位点: pd.DataFrame) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(14, 13), sharex=False)
    for ax, 类型 in zip(axes, 类型顺序):
        数据 = 位点.loc[位点["GFP类型"].eq(类型) & 位点["出现次数"].ge(30)].sort_values("常规位点_1based")
        颜色 = np.where(数据["平均校正后残差"].ge(0), "#2A9D8F", "#D1495B")
        ax.bar(数据["常规位点_1based"], 数据["平均校正后残差"], width=1.0, color=颜色)
        ax.axhline(0, color="#333333", lw=0.8)
        ax.set_title(f"{类型}：出现次数 >= 30 的位点")
        ax.set_ylabel("校正残差")
        ax.grid(axis="y", alpha=0.2)
    axes[-1].set_xlabel("常规位点（1-based）")
    fig.suptitle("各 GFP 类型位点价值分布：已按突变数量校正", fontsize=14)
    保存图片("各类型位点价值分布图.png")


def 绘制候选年度变化(年度变化: pd.DataFrame) -> None:
    if 年度变化.empty:
        return
    展示 = 年度变化.reindex(年度变化["频率变化_2025减2024"].abs().sort_values(ascending=False).index).head(20)
    fig, ax = plt.subplots(figsize=(10, 7))
    颜色 = np.where(展示["频率变化_2025减2024"].ge(0), "#2A9D8F", "#D1495B")
    ax.barh(展示["数据突变编号_0based"], 展示["频率变化_2025减2024"], color=颜色)
    ax.axvline(0, color="#333333", lw=1)
    ax.invert_yaxis()
    ax.set_title("候选序列突变频率年度变化：绝对变化最大的 20 项")
    ax.set_xlabel("频率变化（2025 - 2024）")
    ax.set_ylabel("突变（数据位点为 0-based）")
    ax.grid(axis="x", alpha=0.25)
    保存图片("候选序列年度突变频率变化图.png")


def 取各类型两端(表: pd.DataFrame, 分数字段: str, 每类数量: int = 20) -> tuple[pd.DataFrame, pd.DataFrame]:
    高 = (
        表.sort_values(["GFP类型", 分数字段], ascending=[True, False])
        .groupby("GFP类型", group_keys=False)
        .head(每类数量)
    )
    低 = (
        表.sort_values(["GFP类型", 分数字段], ascending=[True, True])
        .groupby("GFP类型", group_keys=False)
        .head(每类数量)
    )
    return 高, 低


def 写入Excel(数据表: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(工作簿路径, engine="openpyxl") as writer:
        for 工作表, 数据 in 数据表.items():
            数据.to_excel(writer, sheet_name=工作表[:31], index=False)
        for 工作表 in writer.book.worksheets:
            工作表.freeze_panes = "A2"
            工作表.auto_filter.ref = 工作表.dimensions
            for 单元格 in 工作表[1]:
                字体 = copy(单元格.font)
                字体.bold = True
                单元格.font = 字体
            for 列 in 工作表.columns:
                最大宽度 = min(max(len(str(x.value)) if x.value is not None else 0 for x in 列[:300]) + 2, 45)
                工作表.column_dimensions[列[0].column_letter].width = max(最大宽度, 10)


def df_html(数据: pd.DataFrame, 行数: int = 12) -> str:
    return 数据.head(行数).to_html(index=False, border=0, classes="dataframe", float_format=lambda x: f"{x:.4f}")


def 分类型示例_html(数据: pd.DataFrame, 每类行数: int = 5) -> str:
    if 数据.empty:
        return "<p>没有满足严格阈值的记录。</p>"
    示例 = 数据.groupby("GFP类型", group_keys=False, sort=False).head(每类行数)
    return df_html(示例, len(示例))


def 生成报告(
    类型概览: pd.DataFrame,
    数据质量: pd.DataFrame,
    全类型高价值位点: pd.DataFrame,
    全类型低价值位点: pd.DataFrame,
    高价值位点: pd.DataFrame,
    低价值位点: pd.DataFrame,
    候选序列差异: pd.DataFrame,
) -> None:
    图片列表 = [
        ("各类型亮度范围图", "各类型亮度范围图.png"),
        ("各类型亮度分布图", "各类型亮度分布图.png"),
        ("突变数量与亮度关系图", "突变数量与亮度关系图.png"),
        ("各类型位点价值分布图", "各类型位点价值分布图.png"),
        ("候选序列年度突变频率变化图", "候选序列年度突变频率变化图.png"),
    ]
    图片html = "\n".join(
        f"<section><h2>{html.escape(标题)}</h2><img src='图表/{html.escape(文件名)}' alt='{html.escape(标题)}'></section>"
        for 标题, 文件名 in 图片列表
        if (图表目录 / 文件名).exists()
    )
    内容 = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>GFP 综合分析报告</title>
<style>
body {{ font-family: "Microsoft YaHei", Arial, sans-serif; margin: 32px auto; max-width: 1280px; color: #263238; line-height: 1.6; padding: 0 20px; }}
h1, h2 {{ color: #174A5B; }}
h1 {{ border-bottom: 3px solid #2A9D8F; padding-bottom: 10px; }}
h2 {{ margin-top: 34px; border-bottom: 1px solid #ccd9dc; }}
.note {{ background: #f1f7f7; border-left: 4px solid #2A9D8F; padding: 12px 16px; }}
.warn {{ background: #fff7ed; border-left: 4px solid #E76F51; padding: 12px 16px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; overflow-wrap: anywhere; }}
th, td {{ border: 1px solid #d8e0e2; padding: 6px 8px; text-align: left; }}
th {{ background: #e9f3f3; position: sticky; top: 0; }}
img {{ max-width: 100%; border: 1px solid #d8e0e2; }}
code {{ background: #eef2f2; padding: 2px 4px; }}
</style>
</head>
<body>
<h1>GFP 综合分析报告</h1>
<p class="note">完整明细见 <code>GFP综合分析结果.xlsx</code>。报告中的高价值和低价值位点基于关联分析，用于候选筛选和实验优先级排序，不代表已经证明的因果关系。</p>
<h2>分析口径</h2>
<ul>
<li><code>brightness</code> 中的突变位置采用 0-based 编号，例如 <code>A109D</code> 对应参考序列第 110 位。</li>
<li>位点和替换关联价值使用按 GFP 类型及突变数量校正后的亮度残差，降低多点突变数量造成的混杂。</li>
<li>全类型分析先在每种 GFP 内计算分数，再跨类型等权汇总；分类型分析则只在单一 GFP 类型内部排序。</li>
<li>候选序列与 <code>avGFP</code> 参考序列比较，同时保留相对 2024 年首条序列的临时参考差异。</li>
</ul>
<h2>类型概览</h2>
{df_html(类型概览, 20)}
<h2>数据质量摘要</h2>
{df_html(数据质量, 30)}
<h2>全类型高价值位点示例</h2>
<p class="warn">全类型位点按各参考序列的常规编号合并，用于跨 GFP 类型趋势筛选；不同 GFP 的同一编号不保证是严格结构同源位点。</p>
{df_html(全类型高价值位点, 20)}
<h2>全类型低价值位点示例</h2>
{df_html(全类型低价值位点, 20)}
<h2>高价值位点示例</h2>
<p>以下表格为每一种 GFP 类型内部单独计算的结果。</p>
{分类型示例_html(高价值位点)}
<h2>低价值位点示例</h2>
{分类型示例_html(低价值位点)}
<h2>候选序列差异</h2>
{df_html(候选序列差异.drop(columns=["序列"]), 20)}
{图片html}
<h2>补充图表</h2>
<p>更多图表位于 <code>图表</code> 目录。</p>
</body>
</html>
"""
    报告路径.write_text(内容, encoding="utf-8")


def main() -> None:
    设置绘图样式()
    if 输出目录.exists():
        shutil.rmtree(输出目录)
    图表目录.mkdir(parents=True)

    亮度, 候选, 参考序列, 展开 = 构造基础数据()
    数据质量 = 生成数据质量表(亮度, 候选, 参考序列, 展开)
    类型概览, 亮度分位数, 突变数量统计 = 生成类型概览(亮度)
    替换关联, 位点关联, 替换类型 = 生成替换与位点关联表(展开)
    全类型位点价值, 全类型替换类型价值 = 生成全类型价值表(位点关联, 替换类型)
    共现组合 = 生成共现组合表(亮度)
    候选序列差异, 候选突变明细, 候选年度频率, 候选年度变化 = 生成候选序列分析(候选, 参考序列)

    可靠替换 = 替换关联.loc[替换关联["出现次数"].ge(10)]
    高价值替换, _ = 取各类型两端(
        可靠替换.loc[可靠替换["价值分类"].eq("高价值")], "平均校正后残差"
    )
    _, 低价值替换 = 取各类型两端(
        可靠替换.loc[可靠替换["价值分类"].eq("低价值")], "平均校正后残差"
    )
    可靠位点 = 位点关联.loc[位点关联["出现次数"].ge(30)]
    高价值位点, _ = 取各类型两端(
        可靠位点.loc[可靠位点["价值分类"].eq("高价值")], "平均校正后残差"
    )
    _, 低价值位点 = 取各类型两端(
        可靠位点.loc[可靠位点["价值分类"].eq("低价值")], "平均校正后残差"
    )
    全类型高价值位点 = 全类型位点价值.loc[
        全类型位点价值["跨类型价值分类"].eq("跨类型高价值")
    ].sort_values("平均校正后残差_类型等权", ascending=False)
    全类型低价值位点 = 全类型位点价值.loc[
        全类型位点价值["跨类型价值分类"].eq("跨类型低价值")
    ].sort_values("平均校正后残差_类型等权", ascending=True)
    全类型分歧位点 = 全类型位点价值.loc[
        全类型位点价值["跨类型价值分类"].eq("跨类型分歧")
    ].sort_values("平均校正后残差_类型等权", ascending=False)
    全类型高价值替换类型 = 全类型替换类型价值.loc[
        全类型替换类型价值["跨类型价值分类"].eq("跨类型高价值")
    ].sort_values("平均校正后残差_类型等权", ascending=False)
    全类型低价值替换类型 = 全类型替换类型价值.loc[
        全类型替换类型价值["跨类型价值分类"].eq("跨类型低价值")
    ].sort_values("平均校正后残差_类型等权", ascending=True)
    高频共现组合 = 共现组合.sort_values(["GFP类型", "共同出现次数"], ascending=[True, False])
    价值分类统计 = pd.concat(
        [
            位点关联.groupby(["GFP类型", "价值分类"]).size().reset_index(name="记录数").assign(分析维度="位点关联"),
            替换关联.groupby(["GFP类型", "价值分类"]).size().reset_index(name="记录数").assign(分析维度="替换关联"),
            全类型位点价值.groupby("跨类型价值分类").size().reset_index(name="记录数").assign(分析维度="全类型位点", GFP类型="全类型").rename(columns={"跨类型价值分类": "价值分类"}),
            全类型替换类型价值.groupby("跨类型价值分类").size().reset_index(name="记录数").assign(分析维度="全类型替换类型", GFP类型="全类型").rename(columns={"跨类型价值分类": "价值分类"}),
        ],
        ignore_index=True,
    )[["分析维度", "GFP类型", "价值分类", "记录数"]]

    绘制类型亮度范围(亮度)
    绘制类型亮度分布(亮度)
    绘制突变数量关系(突变数量统计)
    绘制位点价值(位点关联)
    绘制候选年度变化(候选年度变化)

    说明 = pd.DataFrame(
        [
            {"项目": "输入文件", "内容": str(输入文件)},
            {"项目": "分析时间口径", "内容": "由脚本运行时生成；不依赖在线数据"},
            {"项目": "位点编号", "内容": "原始 brightness 数据为 0-based；结果同时提供常规 1-based 位点"},
            {"项目": "高低价值位点", "内容": "按同 GFP 类型、同突变数量校正后的平均亮度残差排序；位点至少出现 30 次"},
            {"项目": "全类型高低价值位点", "内容": "每种 GFP 内先计算位点校正残差，再按常规位点编号跨类型等权汇总；同编号不保证严格结构同源"},
            {"项目": "高低价值替换", "内容": "按同 GFP 类型、同突变数量校正后的平均亮度残差排序；具体替换至少出现 10 次"},
            {"项目": "全类型替换类型", "内容": "按 A>G 等替换类型跨 GFP 类型等权汇总；分类要求覆盖至少 2 种 GFP 且总出现次数不少于 20"},
            {"项目": "共现组合", "内容": "仅表示共同出现序列的关联特征，不应解释为已证明的协同效应"},
            {"项目": "候选序列", "内容": "只做序列差异和年度频率分析"},
        ]
    )
    Excel数据表 = {
        "说明": 说明,
        "数据质量": 数据质量,
        "类型概览": 类型概览,
        "亮度分位数": 亮度分位数,
        "突变数量统计": 突变数量统计,
        "价值分类统计": 价值分类统计,
        "全类型位点价值": 全类型位点价值,
        "全类型高价值位点": 全类型高价值位点,
        "全类型低价值位点": 全类型低价值位点,
        "全类型分歧位点": 全类型分歧位点,
        "位点价值总表": 位点关联,
        "高价值位点": 高价值位点,
        "低价值位点": 低价值位点,
        "全类型替换类型价值": 全类型替换类型价值,
        "全类型高价值替换类型": 全类型高价值替换类型,
        "全类型低价值替换类型": 全类型低价值替换类型,
        "替换关联价值": 替换关联,
        "高价值替换": 高价值替换,
        "低价值替换": 低价值替换,
        "替换类型概览": 替换类型,
        "高频共现组合": 高频共现组合,
        "共现组合总表": 共现组合,
        "候选序列差异": 候选序列差异,
        "候选突变明细": 候选突变明细,
        "候选年度频率": 候选年度频率,
        "候选年度变化": 候选年度变化,
    }
    写入Excel(Excel数据表)
    生成报告(
        类型概览,
        数据质量,
        全类型高价值位点,
        全类型低价值位点,
        高价值位点,
        低价值位点,
        候选序列差异,
    )

    print(f"已生成工作簿: {工作簿路径}")
    print(f"已生成报告: {报告路径}")
    print(f"图表数量: {len(list(图表目录.glob('*.png')))}")
    print("类型概览:")
    print(类型概览.to_string(index=False))


if __name__ == "__main__":
    main()
