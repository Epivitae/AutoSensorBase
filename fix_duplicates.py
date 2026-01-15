import json
import os

PROCESSED_FILE = "processed_probes.json"

def clean_duplicates():
    print(f"🧹 开始清理 {PROCESSED_FILE} 中的重复项...")

    # 1. 检查文件是否存在
    if not os.path.exists(PROCESSED_FILE):
        print(f"❌ 文件不存在: {PROCESSED_FILE}")
        return

    # 2. 读取数据
    with open(PROCESSED_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    initial_count = len(data)
    print(f"📚 当前总条数: {initial_count}")
    
    unique_data = []
    seen_dois = set()
    duplicates_count = 0
    
    # 3. 遍历并去重
    for item in data:
        # 获取 DOI (兼容 'doi' 或 'DOI' 键名)
        raw_doi = item.get('doi') or item.get('DOI')
        
        # 如果没有 DOI，为了安全起见，我们保留它（或者你可以选择根据标题去重）
        if not raw_doi:
            unique_data.append(item)
            continue
            
        # === 核心：DOI 标准化 ===
        # 1. 转小写
        # 2. 去除首尾空格
        # 3. 去除 URL 前缀 (兼容 https 和 http)
        clean_doi = raw_doi.lower().strip()
        clean_doi = clean_doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
        
        if clean_doi in seen_dois:
            # 发现重复！跳过不存
            duplicates_count += 1
            # 可选：打印出被删除的重复项名字，方便确认
            print(f"   🗑️  删除重复项: {clean_doi} ({item.get('probe_name', 'Unknown Probe')})")
            continue
        else:
            # 第一次见到，加入集合并保存
            seen_dois.add(clean_doi)
            unique_data.append(item)
            
    # 4. 保存回文件
    if duplicates_count > 0:
        with open(PROCESSED_FILE, "w", encoding="utf-8") as f:
            json.dump(unique_data, f, indent=4, ensure_ascii=False)
        print("=" * 40)
        print(f"✨ 清理完成！")
        print(f"🔻 移除了: {duplicates_count} 条重复记录")
        print(f"📚 剩余: {len(unique_data)} 条有效记录")
    else:
        print("✅ 文件很干净，没有发现重复项。")

if __name__ == "__main__":
    clean_duplicates()