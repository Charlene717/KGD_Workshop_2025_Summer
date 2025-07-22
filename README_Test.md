# 蟹足腫與皮膚遺傳團隊工作坊  
# Workshop of Keloid and Genodermatosis Study  

## 高通量定序分析技術實務應用：單細胞 RNA 定序與空間轉錄組  
### Practical Applications of High‑Throughput Sequencing Technologies: Single‑Cell RNA Sequencing and Spatial Transcriptomics  

- **時間 | Date:** 2025/7/18 (五) 21:00 – 23:00 | July 18 2025 (Fri) 21:00–23:00  
- **地點 | Venue:** 線上會議 [Google Meet](https://meet.google.com/agz-dcee-kfa)

| 時間 Time | 報告者 Presenter | 主題 Topic |
|-----------|-----------------|-----------|
| 21:00–21:20 | 黃道揚<br>Daw‑Yang Hwang | 單細胞定序與空間轉錄組平台介紹<br>**Introduction to Single‑Cell RNA Sequencing and Spatial Transcriptomics Platforms** |
| 21:20–21:30 | 巫政霖<br>Cheng‑Lin Wu | 檢體處理流程簡介<br>**Overview of Sample Processing Workflow** |
| 21:30–21:50 | 劉宗霖 Tsung‑Lin Liu<br>林鉎嵃 Sern‑Yan Lim | 單細胞 RNA 定序資料之生物資訊分析流程<br>**Bioinformatics Workflow for Single‑Cell RNA Sequencing** |
| 21:50–22:10 | 蘇柏嵐<br>Po‑Lan Su | 空間轉錄組學：實驗技術操作與資料分析概念<br>**Spatial Transcriptomics: Wet‑Lab Procedures and Analytical Concepts** |
| 22:10–22:25 | Joanne Jerenice J. Añonuevo | 蟹足腫案例實作經驗分享<br>**Keloid Case Study: Practical Implementation Experience** |
| 22:25–22:40 | 許念芸<br>Nien‑Yun Sharon Hsu | 蟹足腫跨物種整合實作經驗分享<br>**Cross‑Species Integration in Keloid: Practical Implementation Experience** |
| 22:40–23:00 | 張嘉容<br>Chia‑Jung Charlene Chang | 實驗室分析資源導覽與 AI 輔助生物資訊分析技巧<br>**Overview of Dry‑Lab Resources and AI‑Assisted Bioinformatics Techniques** |



---

# scRNA‑seq 分析的實際演示和工作流程講解  
# Practical Demonstration & Workflow Explanation of scRNA‑seq Analysis  

- **時間 | Date:** 2025/8/1 (五) 21:00 – 23:00 | Aug 1 2025 (Fri) 21:00–23:00  
- **地點 | Venue:** 線上會議 [Google Meet](https://meet.google.com/fre-vxzi-bae)  
- **講者 | Speakers:** 蘇柏嵐 Po‑Lan Su、張嘉容 Chia‑Jung Chang (Charlene)

| 時間 Time | 主題 Topic |
|-----------|-----------|
| 21:00–21:10 | 單細胞分析流程總覽及核心概念導讀<br>**Overview of scRNA‑seq Workflow with Key Concepts** |
| 21:10–21:15 | scRNA‑seq 資料下載與預處理策略<br>**Strategies for Downloading & Pre‑processing scRNA‑seq Data** |
| 21:15–21:30 | R 語言基礎導論<br>**Foundations of R Programming** |
| 21:30–21:45 | 品質控制方法、前處理與跨樣本資料整合<br>**Quality Control, Pre‑processing & Cross‑Sample Integration** |
| 21:45–22:00 | 細胞類型註解與分群<br>**Cell‑Type Annotation & Clustering** |
| 22:00–22:15 | 差異基因表現與功能富集分析<br>**Differential Gene Expression & Functional Enrichment** |
| 22:15–22:25 | 細胞間通訊分析<br>**Cell‑Cell Communication Analysis** |
| 22:25–22:35 | 細胞軌跡與擬時序分析方法<br>**Trajectory Inference & Pseudotime Analysis** |
| 22:35–22:50 | 結果評估與分析驗證方法<br>**Evaluation of Analytical Results & Validation Techniques** |
| 22:50–23:00 | 如何依研究目標選擇並導入合適的分析工具與方法<br>**Choosing & Implementing Tools Based on Research Objectives** |



---

## 實作事前準備 | Pre‑Workshop Preparation  

1. **電腦與軟體 | Computer & Software**  
   - 安裝 **R ≥ 4.1.3**（建議 4.1.3）與 **RStudio**  
   - 下載網址: <https://posit.co/download/rstudio-desktop/>

2. **範例資料與程式碼 | Example Data & Scripts**  
   - 範例資料: <https://reurl.cc/NYzGM9>  
   - GitHub 倉庫: <https://github.com/KGDLab/KGD_Workshop_2025_Summer>  
   - **套件安裝 (choose one):**  
     - 最新版本：`source("Install_required_packages.R")`  
     - 指定版本：`source("Install_required_packages_KGD_Lab.R")`  
   - 安裝/技術問題請聯絡 Charlene：<mailto:p88071020@gs.ncku.edu.tw>

3. **GitHub**  
   - 建立個人帳號並提供使用者名稱以加入 KGD_Lab  
   - 組織連結: <https://github.com/KGDLab>



---

# 實作作業說明 | Practical Assignment Instructions  

### 重要時程  
- **書面報告繳交期限:** 工作坊結束後 **1 個月內**  
- **口頭報告:** 2025/8/30 (六) 10:00–16:00 (TBD)  
  - **地點:** 成大醫學院 6F 82‑0624  

### 作業內容  
下載「**線上公開的皮膚相關 scRNA‑seq 資料**」並完成完整分析。  
共 **9 + 1** 項檢核 (Checkpoints) — **#0 與 #5–8 為必答**，其餘可視分析需求選答。  
請依序回答下列問題並附上 **圖表、統計量或文獻** 作為佐證。  

| # | 檢核主題 Topic | 檢查點與關鍵問題 Checkpoints & Key Questions |
|---|---------------|---------------------------------------------|
| 0 | **資料來源與預處理記錄**<br>Data Source & Pre‑processing | • 詳細紀錄資料來源、樣本說明、分析流程、R/套件版本、參數設定 |
| 1 | **品質控制 (QC) 設定**<br>Quality Control | • 閾值與過濾策略是否合理？<br>• 是否過度過濾而排除關鍵細胞？ |
| 2 | **主成分數量 (PCs) 選擇**<br>Principal Components | • 依據哪些指標/概念決定 PC 數？ |
| 3 | **批次效應校正**<br>Batch‑Effect Correction | • 整合後是否改善 batch effect 並保留生物訊號？<br>• 跨平台/實驗室/物種整合須注意何事？ |
| 4 | **群集解析度調整**<br>Clustering Resolution | • 解析度是否合理？判斷過高/過低的依據？<br>• 是否使用量化指標優化解析度？ |
| 5 | **細胞類型標註** *(必答／Required)* | • 標註是否符合已知生物學？<br>• 是否出現族群變化或新穎細胞？ |
| 6 | **差異基因表現 (DEG)** *(必答／Required)* | • 篩選門檻與統計方法是否恰當？為什麼？<br>• DEG 與疾病機轉/治療靶點/細胞功能之關聯？ |
| 7 | **功能富集 & 細胞通訊** *(必答／Required)* | • 富集/通訊分析提出哪些機制假說？<br>• 與臨床或病理的聯繫為何？ |
| 8 | **結果整合** *(必答／Required)* | • 各分析結果是否一致？<br>• 是否能整合為一個生物機制模型？ |
| 9 | **軌跡分析** *(若適用／Optional)* | • 如何界定起點 (root)？<br>• 起點選擇對結果有何影響？ |

### 提交格式 | Submission Format  

| 類型 | 要求 |
|------|------|
| **書面報告** | Word 或 PDF 檔，內含**完整程式碼**、**主要圖表**與**結果說明** |
| **口頭報告** | 10‑分鐘簡報 + 5‑分鐘 Q&A |

> 若對作業有任何疑問，請隨時聯絡講者或助教團隊。  
