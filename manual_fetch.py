import json
import os
from data_fetcher import fetch_papers_by_date_range

RAW_FILE = "raw_papers.json"

def main():
    # ===在此处修改你想抓取的日期范围===
    TARGET_START = "2010/01/01"
    TARGET_END = "2019/12/31"
    # ===============================

    print(f"🚀 手动抓取模式: {TARGET_START} 至 {TARGET_END}")

    # 1. 读取现有的 Raw 数据 (用于去重)
    if os.path.exists(RAW_FILE):
        with open(RAW_FILE, "r", encoding="utf-8") as f:
            raw_data = json.load(f)
    else:
        raw_data = []

    # 建立 DOI 索引
    existing_dois = set()
    for item in raw_data:
        d = item.get('doi')
        if d: existing_dois.add(d.lower().strip().replace("https://doi.org/", ""))
    
    print(f"📚 本地已有 Raw 文献: {len(raw_data)} 篇")

    # 2. 执行抓取
    new_candidates = fetch_papers_by_date_range(TARGET_START, TARGET_END)
    
    # 3. 入库 (Deduplication & Staging)
    added_count = 0
    for p in new_candidates:
        clean_doi = p['doi'].lower().strip().replace("https://doi.org/", "")
        
        # 只有当 DOI 不存在时才添加
        if clean_doi not in existing_dois:
            # 初始化状态
            p['ai_analyzed'] = False 
            p['is_probe'] = False
            
            raw_data.append(p)
            existing_dois.add(clean_doi)
            added_count += 1
    
    # 4. 保存
    if added_count > 0:
        with open(RAW_FILE, "w", encoding="utf-8") as f:
            json.dump(raw_data, f, indent=4, ensure_ascii=False)
        print(f"\n✅ 成功添加 {added_count} 篇新文献到 {RAW_FILE}！")
        print("💡 下一步: 请运行 python update_daily.py 开始 AI 分析。")
    else:
        print("\n💤 没有发现新文献 (所有抓取到的都在库里了)。")

if __name__ == "__main__":
    main()