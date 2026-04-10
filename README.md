# 🧬 SmartPrimer | 引物设计综合平台 

![Python](https://img.shields.io/badge/Python-3.8+-0E1117.svg?style=flat-square&logo=python)
![Streamlit](https://img.shields.io/badge/Streamlit-1.20+-FF4B4B.svg?style=flat-square&logo=streamlit)
![Biopython](https://img.shields.io/badge/Engine-Biopython_|_Primer3-00E5FF.svg?style=flat-square)
![UI](https://img.shields.io/badge/UI_Theme-Dark_Mode-161B22.svg?style=flat-square)

**SmartPrimer** 是一款专为现代分子生物学和合成生物学打造的高自动化引物设计 Web 工具。抛弃繁琐的手动设计与格式转换，本平台为您提供从序列读取到批量订购的无缝体验。

🌐 **[在线体验 SmartPrimer]** *https://smartprimer-zdh.streamlit.app/*

---

## ✨ 核心硬核特性 (Key Features)

### 🔬 1. 引擎级多模式兼容
* **Tai Chi Assembly (太极组装):** 硬核支持多达 16 个片段的复杂环化组装计算。
* **Gibson Assembly / Overlap PCR:** 经典无缝克隆与线性拼接，自动计算同源臂与 Tm 值。
* **高精度 qPCR:** 严格排查二聚体与发夹结构，自动控制扩增子长度与引物 Tm 差。

### 🧬 2. SnapGene 格式“黑科技”解析
* **原生支持 `.dna` 文件:** 引入 `Biopython` 引擎，支持直接拖拽 SnapGene 专有格式。
* **智能元件拆解:** 自动剥离并提取 `.dna` 文件中标记的 `Promoter`, `CDS`, `Terminator` 等功能元件序列，直接化为可用片段，告别手动复制粘贴！
* **实时长度追踪:** 序列框底层注入实时计数器，精确显示碱基对 (bp) 长度，保障数据完整性。

### 🛒 3. 工业级“引物购物车”系统
* **全局状态接力:** 切换不同质粒或设计模式，计算结果自动累加至底部购物车，绝不丢失心血。
* **所见即所得编辑:** 基于强交互数据表 (`st.data_editor`)，支持双击直接修改引物序列与备注，精准勾选删除废弃引物，容错率拉满。
* **SnapGene 完美兼容导出:** 一键合并导出标准三列（`名称\t序列\t备注`）的 `.txt` 纯文本，无缝导入 SnapGene 或直接在金唯智等合成公司粘贴下单。

### 🛡️ 4. 实验溯源与防雷机制
* **模板级命名追踪:** 独创 `质粒名-片段数-序号-长度k-PCR模板-方向` 命名规范，彻底杜绝实验中加错模板的致命失误。
* **内切酶雷达:** 自动扫描引物序列，高亮预警 EcoRI, BsaI 等 10 种克隆克星限制性内切酶位点。
  
---

## ⚖️ 声明与致谢 (Acknowledgments & Disclaimer)

> **®️ 商标声明 (Trademark Notice):** > **SnapGene®** 是 Dotmatics 公司的注册商标。本平台 (SmartPrimer) 是一个独立开发的开源辅助工具，与 SnapGene 官方没有任何附属、赞助或合作关系。对 `.dna` 格式的解析支持仅用于科研数据的互操作性 (Interoperability)。
> 
> **⚙️ 引擎致谢 (Engine Acknowledgments):** > 本平台的引物热力学计算核心由权威开源项目 **[Primer3](https://primer3.ut.ee/)** 驱动；底层序列与格式解析引擎由 **[Biopython](https://biopython.org/)** 驱动。在此向伟大的生物信息学开源社区致敬。
> 
> **⚠️ 免责条款 (Disclaimer):** > 本平台完全免费且开源，旨在辅助生命科学研究并提升工作效率。引物设计结果基于热力学理论模型预测，开发者不对下游生物学实验的最终成败承担任何法律或经济责任。建议在实际提交合成订单前，结合具体的实验场景（如高 GC 区域、特殊聚合酶特性等）进行人工复核。

---

## 🚀 快速本地部署 (Quick Start)

想要在本地实验室服务器上运行这套平台？只需极简 3 步：

```bash
# 1. 克隆本仓库
git clone [https://github.com/Zdh129/SmartPrimer.git](https://github.com/Zdh129/SmartPrimer.git)
cd SmartPrimer

# 2. 安装核心生物信息学与前端依赖
pip install streamlit primer3-py pandas biopython

# 3. 启动应用引擎
streamlit run app.py
