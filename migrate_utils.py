import json
import os

def migrate():
    if not os.path.exists("processed_probes.json"):
        print("No processed data found.")
        return

    with open("processed_probes.json", "r", encoding="utf-8") as f:
        processed = json.load(f)

    # 构造 Raw 数据
    raw_list = []
    for item in processed:
        raw_entry = {
            "title": item.get("title") or item.get("Title"),
            "abstract": item.get("abstract") or item.get("Abstract"),
            "doi": item.get("doi") or item.get("DOI"),
            "journal": item.get("journal") or item.get("Journal"),
            "date": item.get("date"),
            "ai_analyzed": True,  # 标记为已处理
            "is_probe": True      # 标记为是探针
        }
        raw_list.append(raw_entry)

    with open("raw_papers.json", "w", encoding="utf-8") as f:
        json.dump(raw_list, f, indent=4, ensure_ascii=False)
    
    print(f"✅ Migrated {len(raw_list)} entries to raw_papers.json")

if __name__ == "__main__":
    migrate()