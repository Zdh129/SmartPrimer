# 🧬 SmartPrimer | 智能核酸引物设计综合平台

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.20+-FF4B4B.svg)
![Primer3](https://img.shields.io/badge/Engine-Primer3-brightgreen.svg)
![License](https://img.shields.io/badge/License-MIT-purple.svg)

**SmartPrimer** 是一款专为现代分子生物学和合成生物学打造的高自动化引物设计 Web 工具。基于 Python 和业界公认的 `primer3` 核心算法，为生命科学研究者提供从常规克隆到复杂多片段无缝组装的“一站式”解决方案。

🌐 **[在线体验 SmartPrimer]** (https://smartprimer-zdh.streamlit.app/)

---

## ✨ 核心功能 (Key Features)

* **🧩 全场景多模式兼容**
    * **Overlap PCR:** 经典线性片段重叠延伸拼接。
    * **Gibson Assembly:** 多片段闭环同源重组，智能计算载体端引物。
    * **Tai Chi Assembly (太极组装):** 硬核支持多达 16 个片段的复杂环化组装。
    * **qPCR 设计:** 严苛的 6 重过滤算法，排除异源二聚体，锁定极短扩增子。
* **📂 智能序列解析引擎**
    * 支持拖拽上传 `.fasta` 或 `.txt` 文件。
    * 自动识别多条序列并无缝填入工作台，彻底告别繁琐的手动复制粘贴。
* **🏷️ 工业级标准化命名**
    * 内置高通量实验室命名规范：`质粒名-总片段数-序号-长度k-方向` (例：`pNew-4-2-1.3k-F`)。
    * 让引物分装、核对与下游实验一目了然。
* **🛡️ 限制性内切酶预警**
    * 智能扫描引物序列，自动排查是否意外引入 EcoRI, BamHI 等 10 种分子克隆常用限制性内切酶位点。
* **📥 一键生成订购单**
    * 自动换算序列长度、计算精准 Tm 值，一键导出无乱码的 CSV 格式订购单，直接发送给合成公司。

---

## 🚀 快速开始 (Quick Start)

如果你希望在本地计算机上运行本平台，请按照以下步骤操作：

### 1. 环境准备
请确保你的电脑上已安装 Python 3.8 或更高版本。建议使用虚拟环境。

### 2. 克隆仓库
```bash
git clone [https://github.com/Zdh129/SmartPrimer.git](https://github.com/Zdh129/SmartPrimer.git)
