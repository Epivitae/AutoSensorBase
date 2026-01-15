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

def normalize_doi(doi_string):
    """统一清洗 DOI 的辅助函数"""
    if not doi_string: return ""
    return doi_string.lower().strip().replace("https://doi.org/", "").replace("http://doi.org/", "")

def main():
    print("🚀 [ETL] Starting Daily Pipeline...")

    # ================= 1. 抓取阶段 =================
    raw_data = load_json(RAW_FILE)
    
    # 建立 Raw 索引
    existing_raw_dois = set()
    for item in raw_data:
        existing_raw_dois.add(normalize_doi(item.get('doi')))
    
    print("🌍 Fetching from PubMed (5 days)...")
    candidates = fetch_broad_probe_papers(days_back=5)
    
    new_raw_count = 0
    for p in candidates:
        clean_doi = normalize_doi(p['doi'])
        if clean_doi not in existing_raw_dois:
            p['ai_analyzed'] = False 
            p['is_probe'] = False
            raw_data.append(p)
            existing_raw_dois.add(clean_doi)
            new_raw_count += 1
            
    if new_raw_count > 0:
        print(f"📥 Staged {new_raw_count} new papers to {RAW_FILE}")
        save_json(RAW_FILE, raw_data)
    else:
        print("💤 No new raw papers found.")

    # ================= 2. 分析阶段 =================
    pending = [p for p in raw_data if not p.get('ai_analyzed')]
    print(f"⏳ Pending Analysis Queue: {len(pending)} papers")
    
    if not pending:
        print("✅ All caught up. Workflow finished.")
        return

    processed_data = load_json(PROCESSED_FILE)
    
    # 建立 Processed 索引 (用于快速查找是否存在)
    processed_dois_set = set()
    for item in processed_data:
        processed_dois_set.add(normalize_doi(item.get('doi')))
    
    # 批处理大小
    BATCH_SIZE = 800 
    batch = pending[:BATCH_SIZE]
    
    analyzed_count = 0
    new_probe_count = 0
    updated_probe_count = 0

    for paper in batch:
        print(f"🤖 Analyzing: {paper['title'][:50]}...")
        try:
            result = analyze_one_paper(paper)
            
            # 无论结果如何，Raw状态都更新为已分析
            paper['ai_analyzed'] = True
            
            if result and result.get('is_new'):
                # 准备最终数据对象
                final_entry = {**paper, **result}
                final_entry.pop('ai_analyzed', None) # 清理掉内部标记
                final_entry.pop('is_probe', None)
                
                # 获取当前 DOI
                current_doi = normalize_doi(paper.get('doi'))
                
                # === 核心修改：覆盖逻辑 ===
                if current_doi in processed_dois_set:
                    print(f"   🔄 UPDATE EXISTING: {result.get('probe_name')} (Overwriting old data)")
                    
                    # 遍历列表找到那个旧的位置，然后替换它
                    for idx, item in enumerate(processed_data):
                        if normalize_doi(item.get('doi')) == current_doi:
                            processed_data[idx] = final_entry # <--- 覆盖！
                            break
                    updated_probe_count += 1
                else:
                    print(f"   🎉 NEW PROBE: {result.get('probe_name')}")
                    processed_data.append(final_entry) # <--- 新增
                    processed_dois_set.add(current_doi)
                    new_probe_count += 1
                
                paper['is_probe'] = True
            else:
                print("   ❌ Rejected")
                paper['is_probe'] = False
            
            analyzed_count += 1
            time.sleep(1) 
            
        except Exception as e:
            print(f"   ⚠️ Analysis Error: {e}")
            continue

    # ================= 3. 保存阶段 =================
    if analyzed_count > 0:
        save_json(RAW_FILE, raw_data)
        save_json(PROCESSED_FILE, processed_data)
        print(f"💾 Saved updates. (+{new_probe_count} new, 🔄 {updated_probe_count} updated)")

if __name__ == "__main__":
    main()
