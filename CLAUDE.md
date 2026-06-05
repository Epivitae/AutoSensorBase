# AutoSensorBase (ASB) — 项目速览

## 一句话定位
全自动遗传编码荧光探针数据库，每日从 PubMed 抓取文献，用 DeepSeek AI 筛选新探针开发论文，以纯静态页面展示。

## 所属团队
中国科学院 (CAS) 神经科学研究所 (ION) — 嵌合探针团队 (Chimera Sensor Team)
联系人: wangk@ion.ac.cn | 官网: www.cns.ac.cn

---

## 核心文件与职责

| 文件 | 角色 |
|------|------|
| `data_fetcher.py` | PubMed 爬虫 — 用 BioPython Entrez API 检索 GEFS 相关文献，支持按天数/日期范围 |
| `data_analyzer.py` | AI 引擎 — 调 DeepSeek Chat API，判断论文是否开发了新探针，提取 probe_name/target/color/type |
| `update_daily.py` | ETL 主流程 — 串联抓取→去重→AI 分析→写入 JSON，带 Rich 进度条 |
| `index.html` | 前端 — 纯静态单页，侧边栏筛选 + 卡片列表，读取 `processed_probes.json` |
| `manual_fetch.py` | 手动补录工具 |
| `fix_duplicates.py` | 去重修复脚本 |

## 数据文件

| 文件 | 大小 | 说明 |
|------|------|------|
| `raw_papers.json` | ~6.4MB | 原始文献池，每篇带 `ai_analyzed` / `is_probe` 标记 |
| `processed_probes.json` | ~1MB | 结构化探针数据（前端直接消费），只含 AI 判定为探针的条目 |
| `update_meta.json` | 几十字节 | `{"last_update": "YYYY-MM-DD HH:MM"}` 北京时间 |

### 探针数据结构
```json
{
  "title": "...",
  "journal": "...",
  "date": "2026",
  "doi": "https://doi.org/...",
  "abstract": "...",
  "is_new": true,
  "reasoning": "短推理依据",
  "probe_name": "iGlucoSnFR2",
  "target": "Glucose",
  "color": "Green",
  "type": "Metabolite Sensor"
}
```

## 自动化

### 1. 每日定时更新 (`daily_update.yml`)
- 触发: UTC 22:00 (北京时间次日 06:00) + workflow_dispatch
- 流程: checkout → pip install → `python update_daily.py` → 上传 artifact → commit & push
- API Key: `DEEPSEEK_API_KEY` 存于 GitHub Secrets

### 2. Issue 众包入口 (`process_probe_issue.yml`)
- 触发: Issue 创建/编辑
- 流程: 解析 Issue body → 正则抓 DOI/摘要 → `process_issue_pipeline()` → 写入 processed_probes.json → 自动 commit + 关闭 Issue

## 技术栈
- Python 3.9+, BioPython (Entrez), Rich (终端美化)
- DeepSeek Chat API (`deepseek-chat` 模型, json_object 模式)
- GitHub Actions (Cron + Issue 事件)
- 纯 HTML/CSS/JS 前端 (无框架，无构建)

## ETL 流程细节 (`update_daily.py`)
1. 加载 `raw_papers.json`，用 DOI 去重
2. 调 `fetch_broad_probe_papers(days_back=5)` 从 PubMed 抓最近 5 天文献
3. 新文献追加到 raw_papers，标记 `ai_analyzed: false`
4. 取前 600 篇未分析的，逐个调 DeepSeek AI
5. AI 判定 `is_new: true` → 写入 `processed_probes.json`（DOI 重复则覆盖）
6. 每 10 篇保存一次 checkpoint，避免中断丢失进度
7. 最后更新 `update_meta.json` 时间戳

## PubMed 检索策略
```
核心: "Genetically encoded" OR "Fluorescent protein" OR "Bioluminescent" ...
功能: "Sensor" OR "Indicator" OR "Probe" OR "Biosensor" ...
已知家族: GCaMP, GECI, GEVI, iSnFR, GRAB, dLight, FRET, BRET ...
排除: Wastewater, Pollutant, Review
```

## AI Prompt 核心规则
- 只看**开发/工程化/优化新探针**的论文，拒绝纯应用和综述
- 颜色推断：-SnFR/-CaMP 后缀默认绿色，R-/jR- 前缀为红色
- 带 guardrail：如果 reasoning 里含 "used"+"existing" 但 AI 判定 is_new=true，自动纠正为 false

## 补充说明
- 原计划用 Streamlit 部署（README 提到），但实际前端是纯静态 `index.html`
- README 里的 badge 链接用的 `YOUR_USERNAME` 占位符，未替换
- Gemini 相关代码已迁移到 DeepSeek（`data_analyzer.py` 注释里还残留 Gemini 字样）
