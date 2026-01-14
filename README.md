<div align="center">

<img src="https://img.icons8.com/fluency/96/dna-helix.png" alt="logo" width="100" height="100" />

# 🧬 AutoSensorBase (ASB)
### 自动追踪全球最新遗传编码荧光探针数据库
Automated Tracking System for Genetically Encoded Fluorescent Sensors

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/)
![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white)
![Build Status](https://img.shields.io/github/actions/workflow/status/YOUR_USERNAME/AutoSensorBase/daily_update.yml?label=Daily%20Update&logo=github)
![Database Size](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fraw.githubusercontent.com%2FYOUR_USERNAME%2FAutoSensorBase%2Fmain%2Fprocessed_probes.json&query=%24.length&label=Probes&color=success)

![Maintainer](https://img.shields.io/badge/Maintainer-Chimera%20Sensor%20Team-ff69b4?style=flat-square&logo=microgenetics)
![Institution](https://img.shields.io/badge/Institution-CAS%20(ION)-006400?style=flat-square&logo=google-scholar&logoColor=white)
[![Website](https://img.shields.io/badge/Visit-www.cns.ac.cn-blue?style=flat-square)](http://www.cns.ac.cn)
![License](https://img.shields.io/github/license/YOUR_USERNAME/AutoSensorBase?style=flat-square)

[English](./README_EN.md) | [简体中文](./README.md)

</div>

---

## 📖 项目简介 (Introduction)

**AutoSensorBase (ASB)** 是一个全自动化的科研情报系统，旨在为神经科学与合成生物学领域的研究者提供最新的工具开发动态。

该系统利用 **GitHub Actions** 每日定时从 PubMed 抓取最新文献，结合 **LLM (Large Language Model)** 智能识别是否为**新开发的遗传编码荧光探针**（如 GCaMP, GRAB, iSnFR 等），并将结构化数据自动更新至网页端。

## ✨ 核心功能 (Features)

* **🕷️ 自动爬虫**: 每日自动检索 PubMed，追踪 `Genetically Encoded Sensor` 相关文献。
* **🧠 AI 筛选**: 集成大语言模型（GLM-4 / DeepSeek），精准剔除综述与纯应用型文章，只保留工具开发类工作。
* **📊 自动建库**: 提取探针名称、检测底物、荧光颜色、类型等关键参数。
* **🚀 云端展示**: 基于 Streamlit 的交互式网页，支持筛选、搜索与一键导出 CSV。

## 📸 界面预览 (Screenshots)

<div align="center">
  <img src="https://via.placeholder.com/800x400.png?text=AutoSensorBase+Dashboard+Screenshot" alt="Dashboard" width="800" />
</div>

## 🛠️ 技术栈 (Tech Stack)

* **Backend**: Python, `BioPython` (Entrez API)
* **AI Engine**: Zhipu AI (GLM-4) / SiliconFlow (DeepSeek-V3)
* **Automation**: GitHub Actions (Cron Job)
* **Frontend**: Streamlit Cloud
* **Data Storage**: JSON (Git-Scraping)

## 🚀 快速开始 (Quick Start)

如果你想在本地运行本项目：

1.  **克隆仓库**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/AutoSensorBase.git](https://github.com/YOUR_USERNAME/AutoSensorBase.git)
    cd AutoSensorBase
    ```

2.  **安装依赖**
    ```bash
    pip install -r requirements.txt
    ```

3.  **配置 API Key**
    在本地环境变量中设置你的 LLM API Key (例如智谱或 DeepSeek)。

4.  **运行网页**
    ```bash
    streamlit run app.py
    ```

## 👥 维护团队 (Maintainers)

本项目由 **中国科学院 (CAS) 嵌合探针团队 (Chimera Sensor Team)** 开发与维护。

我们致力于开发新型遗传编码生物传感器，以解析复杂的神经信号与代谢网络。

* **所属机构**: Institute of Neuroscience (ION), CAS
* **官方网站**: [www.cns.ac.cn](http://www.cns.ac.cn)
* **联系我们**: wangk@ion.ac.cn

## 📄 许可证 (License)

Distributed under the MIT License. See `LICENSE` for more information.

<div align="center">
  <sub>Built with ❤️ by the Scientific Community.</sub>
</div>