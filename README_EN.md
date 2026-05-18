<div align="center">

<img src="/image/pulumi.svg" alt="logo" width="100" height="100" />

# AutoSensorBase (ASB)
### Automated Tracking System for Genetically Encoded Fluorescent Sensors






![Maintainer](https://img.shields.io/badge/Maintainer-Chimera%20Sensor%20Team-ff69b4?style=flat-square&logo=microgenetics)
![Institution](https://img.shields.io/badge/Institution-CAS%20(ION)-006400?style=flat-square&logo=google-scholar&logoColor=white)
[![Website](https://img.shields.io/badge/Visit-www.cns.ac.cn-blue?style=flat-square)](http://www.cns.ac.cn)


[English](./README_EN.md) | [简体中文](./README.md)

</div>

---

## 📖 Introduction

**AutoSensorBase (ASB)** is a fully automated scientific intelligence system designed to provide researchers in neuroscience and synthetic biology with the latest tool development updates.

The system uses **GitHub Actions** to crawl the latest literature from PubMed daily, combined with **LLM (Large Language Model)** to intelligently identify newly developed genetically encoded fluorescent sensors (such as GCaMP, GRAB, iSnFR, etc.), and automatically updates the structured data to the web interface.

## ✨ Core Features

* **🕷️ Automated Crawler**: Daily automatic retrieval of PubMed to track literature related to `Genetically Encoded Sensor`.
* **🧠 AI Filtering**: Integrated with large language models (GLM-4 / DeepSeek) to accurately exclude reviews and pure application articles, retaining only tool development work.
* **📊 Automated Database Construction**: Extracts key parameters such as probe name, detection substrate, fluorescence color, and type.
* **🚀 Cloud Display**: Interactive web interface based on Streamlit, supporting filtering, searching, and one-click CSV export.

## 📸 Screenshots

<div align="center">
  <img src="/image/snap.png" alt="Dashboard" width="600" />
</div>

## 🛠️ Tech Stack

* **Backend**: Python, `BioPython` (Entrez API)
* **AI Engine**: Zhipu AI (GLM-4) / SiliconFlow (DeepSeek-V3)
* **Automation**: GitHub Actions (Cron Job)
* **Frontend**: Streamlit Cloud
* **Data Storage**: JSON (Git-Scraping)

## 🚀 Quick Start

If you want to run this project locally:

1.  **Clone the repository**
    ```bash
    git clone https://github.com/YOUR_USERNAME/AutoSensorBase.git
    cd AutoSensorBase

