import json
import os
import time
from Bio import Entrez
from data_analyzer import analyze_one_paper

# 配置
Entrez.email = "wangk@ion.ac.cn"
DB_FILE = "processed_probes.json"
TARGET_YEAR_START = "2023/01/01"
TARGET_YEAR_END = "2023/12/31"

def fetch_2023_candidates():
    """专门抓取 2023 年全年的数据"""
    print(f"🌍 正在搜索 PubMed ({TARGET_YEAR_START} - {TARGET_YEAR_END})...")
    
    # === 复用你最新的检索策略 ===
    core_tech = '("Genetically encoded"[Title/Abstract] OR "Fluorescent protein"[Title/Abstract] OR "Fluorescent protein-based"[Title/Abstract] OR "Bioluminescent"[Title/Abstract] OR "Chemogenetic"[Title/Abstract])'
    function_terms = '("Sensor"[Title/Abstract] OR "Indicator"[Title/Abstract] OR "Probe"[Title/Abstract] OR "Reporter"[Title/Abstract] OR "Biosensor"[Title/Abstract])'
    specific_families = '("GCaMP" OR "GECI" OR "GEVI" OR "iSnFR" OR "GRAB" OR "dLight" OR "FRET biosensor" OR "BRET biosensor" OR "cpGFP" OR "CaBLAM")'
    noise_filter = 'NOT ("Wastewater" OR "Pollutant" OR "Review"[Publication Type])'

    # 这里的关键是加入日期范围限定
    date_query = f'("{TARGET_YEAR_START}"[Date - Publication] : "{TARGET_YEAR_END}"[Date - Publication])'
    
    search_query = f'(({core_tech} AND {function_terms}) OR {specific_families}) {noise_filter} AND {date_query}'
    
    try:
        # 搜索 ID (这里 retmax 设大一点，防止漏掉)
        handle = Entrez.esearch(db="pubmed", term=search_query, retmax=2000)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        print(f"📦 PubMed 返回了 {len(id_list)} 篇 2023 年的候选文献 ID。")
        
        if not id_list: return []

        # 批量获取详情
        print("📥 正在下载文献元数据...")
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        papers_data = Entrez.read(handle)
        handle.close()
        
        results = []
        dev_keywords = ["develop", "engineer", "design", "novel", "new", "variant", "optimize"]

        for paper in papers_data['PubmedArticle']:
            article = paper['MedlineCitation']['Article']
            title = article.get('ArticleTitle', 'No Title')
            
            doi = "No DOI"
            for eid in article.get('ELocationID', []):
                if eid.attributes.get('EIdType') == 'doi':
                    doi = str(eid)
            
            # 必须要有 DOI 才有价值
            if doi == "No DOI": continue

            abstract_list = article.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            if not abstract: abstract = "[No Abstract]"

            # 简单的本地关键词初筛 (只保留看起来像开发的)
            text = (title + " " + abstract).lower()
            # if any(kw in text for kw in dev_keywords): # 如果想跑全量，可以注释掉这行筛选
            results.append({
                "title": title,
                "abstract": abstract,
                "doi": f"https://doi.org/{doi}",
                "journal": article.get('Journal', {}).get('Title', 'Unknown'),
                "date": "2023"
            })
            
        return results

    except Exception as e:
        print(f"❌ 搜索出错: {e}")
        return []

def main():
    # 1. 读取现有数据库
    if os.path.exists(DB_FILE):
        with open(DB_FILE, "r", encoding="utf-8") as f:
            current_data = json.load(f)
    else:
        current_data = []
    
    # 建立已存在 DOI 的集合 (用于去重)
    existing_dois = set()
    for item in current_data:
        d = item.get('doi') or item.get('DOI')
        if d: existing_dois.add(d.lower().strip().replace("https://doi.org/", ""))

    print(f"📚 本地已有数据: {len(current_data)} 条")

    # 2. 抓取 2023 数据
    candidates = fetch_2023_candidates()
    
    # 3. 筛选出真正需要 AI 跑的新数据
    to_process = []
    for p in candidates:
        clean_doi = p['doi'].lower().replace("https://doi.org/", "").strip()
        if clean_doi not in existing_dois:
            to_process.append(p)
            
    print(f"\n🔍 经过去重，发现 {len(to_process)} 篇 2023 年的新文献需要 AI 分析。")
    
    if not to_process:
        print("✅ 所有文献都已存在，无需更新。")
        return

    # ⚠️ 确认提示 (防止一次性扣太多钱)
    print("⚠️  注意：这将调用 AI 接口进行分析。")
    confirm = input(f"   输入 'y' 开始分析这 {len(to_process)} 篇文献，输入 'n' 退出: ")
    if confirm.lower() != 'y':
        return

    # 4. 开始处理 (带进度保存)
    print("\n🚀 开始 AI 分析流水线...")
    new_probes_count = 0
    
    for i, paper in enumerate(to_process):
        print(f"   [{i+1}/{len(to_process)}] {paper['title'][:60]}...")
        
        try:
            # AI 分析
            result = analyze_one_paper(paper)
            
            if result and result.get('is_new'):
                print(f"      🎉 发现新探针: {result.get('probe_name')}")
                # 合并数据
                merged = {**paper, **result}
                current_data.append(merged)
                new_probes_count += 1
            else:
                print("      Pass (非新探针开发)")
            
            # 频率限制
            time.sleep(1)

            # === 断点续传机制：每 5 篇保存一次 ===
            if (i + 1) % 5 == 0:
                print("      💾 (自动保存进度...)")
                with open(DB_FILE, "w", encoding="utf-8") as f:
                    json.dump(current_data, f, indent=4, ensure_ascii=False)
                    
        except Exception as e:
            print(f"      ❌ 处理出错: {e}")
            continue

    # 5. 最后再一次保存
    with open(DB_FILE, "w", encoding="utf-8") as f:
        json.dump(current_data, f, indent=4, ensure_ascii=False)
    
    print(f"\n🎉 2023数据回填完成！共新增 {new_probes_count} 个探针。")

if __name__ == "__main__":
    main()