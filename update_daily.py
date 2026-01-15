import json
import os
import time
# ✅ 改动1：引入新的广度优先爬虫函数
from data_fetcher import fetch_broad_probe_papers
from data_analyzer import analyze_one_paper

DB_FILE = "processed_probes.json"

def main():
    print("🚀 Starting Daily Update...")

    # ================= 1. 读取旧数据 =================
    if os.path.exists(DB_FILE):
        try:
            with open(DB_FILE, "r", encoding="utf-8") as f:
                old_data = json.load(f)
        except json.JSONDecodeError:
            print("⚠️ JSON file corrupted, starting fresh.")
            old_data = []
    else:
        old_data = []
    
    # 提取现有 DOI 集合（统一转小写，防止大小写差异导致的重复）
    existing_dois = set()
    for item in old_data:
        d = item.get('doi') or item.get('DOI')
        if d:
            existing_dois.add(d.lower().strip())

    print(f"📚 Existing database has {len(old_data)} entries.")

    # ================= 2. 爬取新数据 =================
    # ✅ 改动2：使用 fetch_broad_probe_papers
    # 建议设置 3-5 天，因为我们现在用的是"录入日期(Entrez Date)"，这个更新非常及时
    print("🌍 Fetching data using BROAD strategy (Entrez Date)...")
    try:
        new_papers = fetch_broad_probe_papers(days_back=5)
    except Exception as e:
        print(f"❌ Critical Error during fetching: {e}")
        return

    print(f"📦 Fetched {len(new_papers)} candidates from PubMed.")

    # ================= 3. 筛选未处理的 =================
    to_process = []
    for p in new_papers:
        # 获取 DOI 并清洗
        current_doi = p.get('doi') or p.get('DOI')
        if not current_doi:
            continue
            
        clean_doi = current_doi.lower().strip()
        
        # 只有当 DOI 不在库中时，才处理
        if clean_doi not in existing_dois:
            to_process.append(p)

    print(f"🔍 Found {len(to_process)} NEW papers to analyze (after deduplication).")

    if not to_process:
        print("💤 No new unique papers found. Exiting.")
        return

    # ================= 4. AI 分析 =================
    # 提示：如果一次更新太多（例如 >20 篇），可能需要考虑分批运行
    added_count = 0
    
    print("🤖 Starting AI Analysis...")
    
    for i, paper in enumerate(to_process):
        title = paper.get('title') or paper.get('Title')
        print(f"   [{i+1}/{len(to_process)}] Analyzing: {title[:50]}...")
        
        # 调用 AI
        result = analyze_one_paper(paper)
        
        # ✅ 改动3：增加延时，防止 API 并发过高报错
        time.sleep(1.0) 
        
        if result and result.get('is_new'):
            probe_name = result.get('probe_name', 'Unknown')
            print(f"      ✅ NEW PROBE DISCOVERED: {probe_name}")
            
            # 合并原始数据和 AI 分析结果
            # 注意：保留 paper 里的 metadata，用 result 覆盖关键字段
            merged_entry = {**paper, **result}
            
            # 确保统一的小写 key 存在 (为了兼容前端)
            if 'Title' in merged_entry and 'title' not in merged_entry:
                merged_entry['title'] = merged_entry['Title']
            if 'Abstract' in merged_entry and 'abstract' not in merged_entry:
                merged_entry['abstract'] = merged_entry['Abstract']
                
            old_data.append(merged_entry) # 加入总表
            added_count += 1
        else:
            # 可能是纯应用文章，或者是提取失败
            print("      ❌ Not a new sensor development.")

    # ================= 5. 保存结果 =================
    if added_count > 0:
        print(f"💾 Saving {added_count} new entries to {DB_FILE}...")
        # 备份一下是个好习惯（可选）
        # shutil.copy(DB_FILE, DB_FILE + ".bak") 
        
        with open(DB_FILE, "w", encoding="utf-8") as f:
            json.dump(old_data, f, indent=4, ensure_ascii=False)
        print("🎉 Update Complete!")
    else:
        print("🏁 Analysis complete, no new probes added to database.")

if __name__ == "__main__":
    main()