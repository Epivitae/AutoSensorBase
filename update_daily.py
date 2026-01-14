import json
import os
import time
from data_fetcher import fetch_recent_papers
from data_analyzer import analyze_one_paper

DB_FILE = "processed_probes.json"

def main():
    print("🚀 Starting Daily Update...")

    # 1. 读取旧数据
    if os.path.exists(DB_FILE):
        with open(DB_FILE, "r", encoding="utf-8") as f:
            old_data = json.load(f)
    else:
        old_data = []
    
    existing_dois = set(item['doi'] for item in old_data)
    print(f"📚 Existing database has {len(old_data)} entries.")

    # 2. 爬取新数据 (过去 5 天，确保不漏)
    new_papers = fetch_recent_papers(days_back=5)
    print(f"🌍 Fetched {len(new_papers)} recent papers from PubMed.")

    # 3. 筛选未处理的
    to_process = [p for p in new_papers if p['doi'] not in existing_dois]
    print(f"🔍 Found {len(to_process)} NEW papers to analyze.")

    if not to_process:
        print("💤 No new papers to process. Exiting.")
        return

    # 4. AI 分析
    added_count = 0
    for paper in to_process:
        print(f"🤖 Analyzing: {paper['title'][:30]}...")
        result = analyze_one_paper(paper)
        
        if result and result.get('is_new'):
            print(f"   ✅ NEW PROBE: {result['probe_name']}")
            merged_entry = {**paper, **result}
            old_data.append(merged_entry) # 加入总表
            added_count += 1
        else:
            print("   - Pass")
        
        time.sleep(0.5) # 稍微慢点

    # 5. 保存更新
    if added_count > 0:
        with open(DB_FILE, "w", encoding="utf-8") as f:
            json.dump(old_data, f, indent=4, ensure_ascii=False)
        print(f"🎉 Database updated! Added {added_count} new probes.")
    else:
        print("💨 Analysis complete, but no new probes found.")

if __name__ == "__main__":
    main()