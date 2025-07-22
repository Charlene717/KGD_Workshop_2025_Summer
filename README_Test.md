
# 蟹足腫與皮膚遺傳團隊工作坊  <br>Workshop of Keloid and Genodermatosis Study<br>



## Workshop 1 — 高通量定序分析技術實務應用：單細胞 RNA 定序與空間轉錄組  <br>Practical Applications of High‑Throughput Sequencing Technologies: Single‑Cell RNA Sequencing & Spatial Transcriptomics<br>

- **時間 Date:** 2025 / 07 / 18 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** 線上會議 Online Meeting <待加入影片撥放清單>

### 📅 議程 Agenda
[📑 **議程表下載 (Workshop_2025_0718_Agenda.pdf)**](./Agenda_and_Guidelines/Workshop_2025_0718_Agenda.pdf)

| 時間 Time | 報告者 Presenter | 主題 Topic |
|-----------|-----------------|------------|
| 21:00‑21:20 | 黃道揚 Daw‑Yang Hwang | 單細胞定序與空間轉錄組平台介紹<br>Introduction to Single‑Cell RNA Sequencing and Spatial Transcriptomics Platforms |
| 21:20‑21:30 | 巫政霖 Cheng‑Lin Wu | 檢體處理流程簡介<br>Overview of Sample Processing Workflow |
| 21:30‑21:50 | 劉宗霖 Tsung‑Lin Liu<br>林鉎嵃 Sern‑Yan Lim | 單細胞 RNA 定序資料之生物資訊分析流程<br>Bioinformatics Workflow for Single‑Cell RNA Sequencing |
| 21:50‑22:10 | 蘇柏嵐 Po‑Lan Su | 空間轉錄組學：實驗技術操作與資料分析概念<br>Spatial Transcriptomics: Wet‑Lab Procedures & Analytical Concepts |
| 22:10‑22:25 | Joanne Jerenice J. Añonuevo | 蟹足腫案例實作經驗分享<br>Keloid Case Study: Practical Implementation Experience |
| 22:25‑22:40 | 許念芸 Nien‑Yun Sharon Hsu | 蟹足腫跨物種整合實作經驗分享<br>Cross‑Species Integration in Keloid: Practical Implementation Experience |
| 22:40‑23:00 | 張嘉容 Chia‑Jung Charlene Chang | 實驗室分析資源導覽與 AI 輔助生物資訊分析技巧<br>Overview of Dry‑Lab Resources & AI‑Assisted Bioinformatics Techniques |

---

## Workshop 2 — scRNA‑seq 分析的實際演示與工作流程講解  ·<br>Practical Demonstration and Workflow Explanation of scRNA‑seq Analysis<br>

- **時間 Date:** 2025 / 08 / 01 (五 Friday) 21:00 – 23:00  
- **地點 Venue:** Online Meeting <待加入影片撥放清單>  
- **講者 Speakers:** 蘇柏嵐 Po‑Lan Su · 張嘉容 Chia‑Jung Chang (Charlene)

### 📅 議程 Agenda
[📑 **議程表下載 (Workshop_2025_0801_Agenda.pdf)**](./Agenda_and_Guidelines/Workshop_2025_0801_Agenda.pdf)

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

---

## 實作事前準備 / Pre‑Workshop Preparation

### 1. 電腦與軟體 / Computer & Software
- 安裝 **R ≥ 4.1.3**（實驗室慣用版本 4.1.3）與 **RStudio**。  
  Download: <https://posit.co/download/rstudio-desktop/>

### 2. 範例資料與程式碼 / Example Data & Scripts
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

---

## 實作作業說明 / Practical Assignment Instructions

- **Written report deadline:** within **1 month** after workshop  
- **Oral presentation:** 2025 / 08 / 30 (Saturday) 10:00 – 16:00 (TBD)  
  - Venue: College of Medicine, NCKU — Room 82‑0624

### Assignment Overview
Participants must download a *public skin‑related scRNA‑seq dataset* and perform a **complete analysis** following the checkpoints below.

### 9 + 1 Checkpoints
Checkpoints 0 & 5 – 8 are **mandatory**; others are optional.

| # | 檢核主題 Topic | 檢查點與關鍵問題 Checkpoints & Key Questions |
|---|---------------|---------------------------------------------|
| 0 | 資料來源與預處理記錄<br>Data Source & Pre‑processing | • 詳細紀錄 dataset、samples、workflow、R & package versions、parameters |
| 1 | 品質管制（QC）設定<br>Quality Control | • 閾值與過濾策略是否合理？<br>• 是否過度過濾排除關鍵細胞？ |
| 2 | 主成分數量（PCs）選擇<br>PC Selection | • 依據哪些指標／概念決定 PC 數？ |
| 3 | 批次效應校正<br>Batch Correction | • 是否減少 batch effect 並保留生物訊號？<br>• 跨平台／物種整合需注意何事？ |
| 4 | 群集解析度調整<br>Clustering Resolution | • 解析度是否合宜？如何判斷過高／過低？<br>• 是否使用量化指標優化？ |
| 5 | **細胞類型標註 (必答)**<br>Cell‑Type Annotation **(Required)** | • 標註是否符合已知生物學？<br>• 有無新穎或顯著變動的細胞族群？ |
| 6 | **差異基因表現 (DEG) (必答)**<br>Differential Expression **(Required)** | • 篩選門檻與統計方法是否恰當？為何？<br>• 這些 DEG 與疾病、靶點、功能的關聯？ |
| 7 | **功能富集 & 細胞通訊 (必答)**<br>Enrichment & Cell‑Cell Communication **(Required)** | • 機制假說為何？<br>• 與臨床／病理如何連結？ |
| 8 | **結果整合 (必答)**<br>Integration of Results **(Required)** | • 各分析結果是否一致？<br>• 能否整合成生物機制模型？ |
| 9 | 軌跡分析（若適用）<br>Trajectory (Optional) | • 如何界定 root？<br>• 起點選擇對結果有何影響？ |

### Submission Format
- **Written report:** Word or PDF — include full code, key plots, explanations  
- **Oral presentation:** *10 min* talk + *5 min* Q&A
