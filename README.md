
<div align="center">

# 蟹足腫與皮膚遺傳團隊工作坊  <br>Workshop of Keloid and Genodermatosis Study<br>
</div>

## 工作坊 1 — 高通量定序分析技術實務應用：單細胞 RNA 定序與空間轉錄組  <br>Workshop 1 — Practical Applications of High‑Throughput Sequencing Technologies: Single‑Cell RNA Sequencing & Spatial Transcriptomics<br>



- **時間 Date:** 2025 / 07 / 18 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** 線上會議 Online Meeting

### 📅 議程 Agenda

| 時間 Time | 報告者 Presenter | 主題 Topic | 
|-----------|-----------------|------------|
| 21:00‑21:20 | 黃道揚 Daw‑Yang Hwang | 單細胞定序與空間轉錄組平台介紹<br>Introduction to Single‑Cell RNA Sequencing and Spatial Transcriptomics Platform| 
| 21:20‑21:30 | 巫政霖 Cheng‑Lin Wu | 檢體處理流程簡介<br>Overview of Sample Processing Workflow | 
| 21:30‑21:50 | 劉宗霖 Tsung‑Lin Liu<br>林鉎嵃 Sern‑Yan Lim | 單細胞 RNA 定序資料之生物資訊分析流程<br>Bioinformatics Workflow for Single‑Cell RNA Sequencing | 
| 21:50‑22:10 | 蘇柏嵐 Po‑Lan Su | 空間轉錄組學：實驗技術操作與資料分析概念<br>Spatial Transcriptomics: Wet‑Lab Procedures & Analytical Concepts | 
| 22:25‑22:40 | 許念芸<br>Nien‑Yun Sharon Hsu<br> | 蟹足腫跨物種整合實作經驗分享<br>Cross‑Species Integration in Keloid: Practical Implementation Experience |
| 22:40‑23:00 | 張嘉容<br>Chia‑Jung Charlene Chang<br> | 實驗室分析資源導覽與 AI 輔助生物資訊分析技巧<br>Overview of Dry‑Lab Resources & AI‑Assisted Bioinformatics Techniques |

<br>

---

## 工作坊 2 — scRNA‑seq 分析的實際演示與工作流程講解  ·<br>Workshop 2 —Practical Demonstration and Workflow Explanation of scRNA‑seq Analysis<br>

- **時間 Date:** 2025 / 08 / 01 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** 線上會議 Online Meeting <待加入影片撥放清單>  
- **講者 Speakers:** 蘇柏嵐 Po‑Lan Su · 張嘉容 Chia‑Jung Chang (Charlene)

### 📅 議程 Agenda
[📑 **議程表下載 (Workshop_2025_Summer_Part2_0801_Agenda.pdf)**](./Agenda_and_Guidelines/Workshop_2025_Summer_Part2_0801_Agenda.pdf)

<div align="center">

| 時間 Time | 主題 Topic |
|-----------|-----------|
| 21:00‑21:10 | 單細胞分析流程總覽及核心概念導讀<br>Overview of scRNA‑seq Workflow with Key Concepts |
| 21:10‑21:15 | scRNA‑seq 資料下載與預處理策略<br>Strategies for Downloading & Pre‑processing scRNA‑seq Data |
| 21:15‑21:30 | R 語言基礎導論<br>Foundations of R Programming |
| 21:30‑21:45 | 品質控制方法、前處理與跨樣本資料整合<br>Quality Control & Pre‑processing & Cross‑Sample Integration |
| 21:45‑22:00 | 細胞類型註解與分群<br>Cell‑Type Annotation & Clustering |
| 22:00‑22:15 | 差異基因表現與功能富集分析<br>Differential Gene Expression & Functional Enrichment |
| 22:15‑22:25 | 細胞間通訊分析<br>Cell‑Cell Communication Analysis |
| 22:25‑22:35 | 細胞軌跡與擬時序分析方法<br>Trajectory Inference & Pseudotime Analysis |
| 22:35‑22:50 | 結果評估與Loupe Browser視覺化互動<br>Results Evaluation & Loupe Browser Interactive Visualization |
| 22:50‑23:00 | 根據研究目標選擇與導入適切的分析工具<br>Selecting & Implementing Tools Based on Research Objectives |

</div>
<br>

## 實作事前準備 Pre‑Workshop Preparation

### 1. 電腦與軟體 Computer & Software  
> **注意:** 單細胞 RNA‑seq 及後續整合／軌跡分析會大量佔用記憶體；若資料規模 > 100 k cells，請優先考慮使用 64 GB RAM 以上或雲端／HPC。  
> **Note:** Single‑cell RNA‑seq, downstream integration, and trajectory analyses are memory‑intensive. For datasets larger than 100 k cells, consider machines with ≥ 64 GB RAM or use cloud/HPC resources.

| 規格項目 Spec | **最低配備 Minimum** | **建議配備 Recommended** |
|--------------|----------------------|--------------------------|
| 作業系統 OS | Windows 10 (64‑bit) / macOS 12 / Ubuntu 20.04 | Windows 11 (64‑bit) / macOS 14 / Ubuntu 22.04 |
| 處理器 CPU | 6 cores (Intel / AMD / Apple Silicon) | 8 – 12 cores, multithreaded / high clock speed |
| 記憶體 RAM | 16 GB | 32 GB (≥ 64 GB for large multi‑sample projects) |
| 儲存空間 Disk | SSD, ≥ 50 GB free | NVMe SSD, ≥ 200 GB free |
| GPU (可選 Optional) | — | CUDA‑capable GPU / Apple M‑series for accelerated visualization & deep learning |
| 網路 Internet | Stable connection ≥ 10 Mbps | ≥ 30 Mbps for faster data & package downloads |
| 額外選項 Extras | — | Docker / Singularity and HPC / Cloud accounts for reproducibility & large‑scale computing |

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
[**實作作業指引檔案 (Workshop_2025_Summer_Part3_Practical_Assignment_Guidelines.pdf)**](./Agenda_and_Guidelines/Workshop_2025_Summer_Part3_Practical_Assignment_Guidelines.pdf)  



&nbsp;
---
## 實驗室網站
[**KGD Lab**](https://twkgd.wordpress.com/)
