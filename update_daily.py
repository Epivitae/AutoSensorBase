import json
import os
import time
from data_fetcher import fetch_broad_probe_papers
from data_analyzer import analyze_one_paper

RAW_FILE = "raw_papers.json"
PROCESSED_FILE = "processed_probes.json"

def load_json(filename):
    if os.path.exists(filename):
        with open(filename, "r", encoding="utf-8") as f:
            return json.load(f)
    return []

def save_json(filename, data):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)

def main():
    print("🚀 [ETL] Starting Daily Pipeline...")

    # ================= 1. 抓取阶段 (Fetch to Raw) =================
    raw_data = load_json(RAW_FILE)
    
    # 建立索引防止重复抓取
    existing_dois = set()
    for item in raw_data:
        d = item.get('doi')
        if d: existing_dois.add(d.lower().strip().replace("https://doi.org/", ""))
    
    # 爬取过去 5 天
    print("🌍 Fetching from PubMed (5 days)...")
    candidates = fetch_broad_probe_papers(days_back=5)
    
    new_raw_count = 0
    for p in candidates:
        clean_doi = p['doi'].lower().strip().replace("https://doi.org/", "")
        if clean_doi not in existing_dois:
            # 初始化状态
            p['ai_analyzed'] = False 
            p['is_probe'] = False
            raw_data.append(p)
            existing_dois.add(clean_doi)
            new_raw_count += 1
            
    if new_raw_count > 0:
        print(f"📥 Staged {new_raw_count} new papers to {RAW_FILE}")
        save_json(RAW_FILE, raw_data)
    else:
        print("💤 No new raw papers found.")

    # ================= 2. 分析阶段 (Process Pending) =================
    # 找出所有未分析的
    pending = [p for p in raw_data if not p.get('ai_analyzed')]
    print(f"⏳ Pending Analysis Queue: {len(pending)} papers")
    
    if not pending:
        print("✅ All caught up. Workflow finished.")
        return

    processed_data = load_json(PROCESSED_FILE)
    
    # 批处理限制 (防止 CI 超时)
    BATCH_SIZE = 1000 
    batch = pending[:BATCH_SIZE]
    
    analyzed_count = 0
    new_probe_count = 0

    for paper in batch:
        print(f"🤖 Analyzing: {paper['title'][:50]}...")
        try:
            result = analyze_one_paper(paper)
            
            # 更新 Raw 状态
            paper['ai_analyzed'] = True
            
            if result and result.get('is_new'):
                print(f"   🎉 NEW PROBE: {result.get('probe_name')}")
                paper['is_probe'] = True
                
                # 存入成品库 (清洗掉内部状态字段)
                final_entry = {**paper, **result}
                final_entry.pop('ai_analyzed', None)
                final_entry.pop('is_probe', None)
                processed_data.append(final_entry)
                new_probe_count += 1
            else:
                print("   ❌ Rejected (Review/App)")
                paper['is_probe'] = False
            
            analyzed_count += 1
            time.sleep(1) # Rate limit
            
        except Exception as e:
            print(f"   ⚠️ Analysis Error: {e}")
            continue

    # ================= 3. 保存阶段 =================
    if analyzed_count > 0:
        save_json(RAW_FILE, raw_data)        # 更新状态
        save_json(PROCESSED_FILE, processed_data)  # 更新成品库
        print(f"💾 Saved updates. (+{new_probe_count} probes)")

if __name__ == "__main__":
    main()