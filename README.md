
<div align="center">

# 蟹足腫與皮膚遺傳團隊工作坊  <br>Workshop of Keloid and Genodermatosis Study<br>
</div>

## 工作坊 1 — 高通量定序分析技術實務應用：單細胞 RNA 定序與空間轉錄組  <br>Workshop 1 — Practical Applications of High‑Throughput Sequencing Technologies: Single‑Cell RNA Sequencing & Spatial Transcriptomics<br>



- **時間 Date:** 2025 / 07 / 18 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** 線上會議 Online Meeting <待加入影片撥放清單>

### 📅 議程 Agenda
[📑 **議程表下載 (Workshop_2025_0718_Agenda.pdf)**](./Agenda_and_Guidelines/Workshop_2025_0718_Agenda.pdf)

| 時間 Time | 報告者 Presenter | 主題 Topic | 影片                                                                   |
|-----------|-----------------|------------|------------|
| 21:00‑21:20 | 黃道揚 Daw‑Yang Hwang | 單細胞定序與空間轉錄組平台介紹<br>Introduction to Single‑Cell RNA Sequencing and Spatial Transcriptomics Platform| [Video]https://www.youtube.com/watch?v=K4p7_n8A1yA&list=PLsjV-b9qanZZR0xkebAP-_jaFGwxtnp5d&index=1Q)                       |
| 21:20‑21:30 | 巫政霖 Cheng‑Lin Wu | 檢體處理流程簡介<br>Overview of Sample Processing Workflow |
| 21:30‑21:50 | 劉宗霖 Tsung‑Lin Liu<br>林鉎嵃 Sern‑Yan Lim | 單細胞 RNA 定序資料之生物資訊分析流程<br>Bioinformatics Workflow for Single‑Cell RNA Sequencing |
| 21:50‑22:10 | 蘇柏嵐 Po‑Lan Su | 空間轉錄組學：實驗技術操作與資料分析概念<br>Spatial Transcriptomics: Wet‑Lab Procedures & Analytical Concepts |
| 22:10‑22:25 | Joanne Jerenice J. Añonuevo | 蟹足腫案例實作經驗分享<br>Keloid Case Study: Practical Implementation Experience |
| 22:25‑22:40 | 許念芸<br>Nien‑Yun Sharon Hsu<br> | 蟹足腫跨物種整合實作經驗分享<br>Cross‑Species Integration in Keloid: Practical Implementation Experience |
| 22:40‑23:00 | 張嘉容<br>Chia‑Jung Charlene Chan<br>g | 實驗室分析資源導覽與 AI 輔助生物資訊分析技巧<br>Overview of Dry‑Lab Resources & AI‑Assisted Bioinformatics Techniques |

<br>

---

## 工作坊 2 — scRNA‑seq 分析的實際演示與工作流程講解  ·<br>Workshop 2 —Practical Demonstration and Workflow Explanation of scRNA‑seq Analysis<br>

- **時間 Date:** 2025 / 08 / 01 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** Online Meeting <待加入影片撥放清單>  
- **講者 Speakers:** 蘇柏嵐 Po‑Lan Su · 張嘉容 Chia‑Jung Chang (Charlene)

### 📅 議程 Agenda
[📑 **議程表下載 (Workshop_2025_0801_Agenda.pdf)**](./Agenda_and_Guidelines/Workshop_2025_0801_Agenda.pdf)

<div align="center">

| 時間 Time | 主題 Topic |
|-----------|-----------|
| 21:00‑21:10 | 單細胞分析流程總覽及核心概念導讀<br>Overview of scRNA‑seq Workflow with Key Concepts |
| 21:10‑21:15 | scRNA‑seq 資料下載與預處理策略<br>Strategies for Downloading & Pre‑processing scRNA‑seq Data |
| 21:15‑21:30 | R 語言基礎導論<br>Foundations of R Programming |
| 21:30‑21:45 | 品質控制方法、前處理與跨樣本資料整合<br>Quality Control, Pre‑processing & Cross‑Sample Integration |
| 21:45‑22:00 | 細胞類型註解與分群<br>Cell‑Type Annotation & Clustering |
| 22:00‑22:15 | 差異基因表現與功能富集分析<br>Differential Gene Expression & Functional Enrichment |
| 22:15‑22:25 | 細胞間通訊分析<br>Cell‑Cell Communication Analysis |
| 22:25‑22:35 | 細胞軌跡與擬時序分析方法<br>Trajectory Inference & Pseudotime Analysis |
| 22:35‑22:50 | 結果評估與分析驗證方法<br>Evaluation & Validation of Analytical Results |
| 22:50‑23:00 | 根據研究目標選擇與導入適切的分析工具<br>Selecting & Implementing Tools Based on Research Objectives |

</div>
<br>

## 實作事前準備 Pre‑Workshop Preparation

### 1. 電腦與軟體 Computer & Software  
> **說明 Note:** 單細胞 RNA-seq 及後續整合／軌跡分析會大量佔用記憶體；若資料規模 > 100k cells，請優先考慮使用 64 GB RAM 以上或雲端／HPC。

| 規格項目 Spec | **最低配備 Minimum** | **建議配備 Recommended** |
|---------------|---------------------|--------------------------|
| 作業系統 OS | Windows 10 (64-bit) / macOS 12 / Ubuntu 20.04 | Windows 11 (64-bit) / macOS 14 / Ubuntu 22.04 |
| 處理器 CPU | 6 cores (Intel / AMD / Apple Silicon) | 8 – 12 cores，多執行緒 / 高時脈 |
| 記憶體 RAM | 16 GB | 32 GB （大量樣本建議 64 GB +） |
| 儲存空間 Disk | SSD，≥ 50 GB free | NVMe SSD，≥ 200 GB free |
| GPU (可選 Optional) | — | CUDA-capable / Apple M-series；加速可視化與深度學習 |
| 網路 Internet | 穩定連線 ≥ 10 Mbps | ≥ 30 Mbps 以便快速下載資料與套件 |
| 額外選項 Extras | — | Docker / Singularity 與 HPC/Cloud 帳號，用於重現及大型計算 |

**English Summary**

- **Minimum:** 6-core CPU, 16 GB RAM, SSD ≥ 50 GB.  
- **Recommended:** 8–12-core CPU, 32 GB RAM (64 GB + for big datasets), NVMe SSD ≥ 200 GB, optional CUDA/Metal-capable GPU for visualization & ML.  
- A stable ≥ 30 Mbps network and access to Docker/Singularity or HPC/Cloud resources are strongly encouraged for large-scale analyses.

---

安裝 **R ≥ 4.1.3**（實驗室慣用版本 4.1.3）與 **RStudio**  
Install **R ≥ 4.1.3** (lab’s preferred 4.1.3) and **RStudio**  
<https://posit.co/download/rstudio-desktop/>


### 2. 範例資料與程式碼 Example Data & Scripts
- Example data: <https://reurl.cc/NYzGM9>  
- GitHub repository (scripts & installers): <https://github.com/KGDLab/KGD_Workshop_2025_Summer>
  ```r
  # Install latest package versions
  source("Install_required_packages.R")

  # OR: install lab‑specific versions
  source("Install_required_packages_KGD_Lab.R")
  ```
- Installation / technical issues → contact **Charlene** <p88071020@gs.ncku.edu.tw>

### 3. GitHub
- Create a personal GitHub account and share your **username** to join **KGD_Lab**.  
  Organization URL: <https://github.com/KGDLab>

&nbsp;
---
## 🛠️ 實作流程導覽 Workshop Practical Workflow
[**實作流程導覽 (Workshop_Practical_Workflow.md)**](./Agenda_and_Guidelines/Workshop_Practical_Workflow.md)  

*本文件逐步示範 scRNA‑seq 分析流程，包含教學腳本與程式碼。*  
*This document walks through the scRNA‑seq analysis workflow with step‑by‑step scripts and code.*

&nbsp;
---
## 📑 實作作業指引 Practical Assignment Guidelines
[**實作作業指引檔案 (Workshop_2025_Practical_Assignment_Guidelines.pdf)**](./Agenda_and_Guidelines/Workshop_2025_Practical_Assignment_Guidelines.pdf)  

*本文件詳述作業目的、繳交格式、9 + 1 項 Checkpoints 及評分標準。*  
*This document outlines the assignment objectives, submission format, 9 + 1 checkpoints, and grading criteria.*


&nbsp;
---
## 實驗室網站
[KGD Lab](https://twkgd.wordpress.com/)
