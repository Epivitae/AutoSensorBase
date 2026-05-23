import time
from Bio import Entrez
from datetime import datetime, timezone, timedelta
import json

# 设置邮箱
Entrez.email = "wangk@ion.ac.cn"


def save_update_timestamp():
    # 设定为北京时间 (UTC+8)
    tz_utc_8 = timezone(timedelta(hours=8))
    current_time = datetime.now(tz_utc_8).strftime("%Y-%m-%d %H:%M")
    
    meta_info = {
        "last_update": current_time
    }
    
    with open("update_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta_info, f, ensure_ascii=False, indent=2)
        
    print(f"✅ Timestamp updated: {current_time}")


def fetch_broad_probe_papers(days_back=5):
    """
    爬取过去 N 天录入 (Entrez Date) 的相关文献
    """
    # 1. 检索词策略
    core_tech = '("Genetically encoded"[Title/Abstract] OR "Fluorescent protein"[Title/Abstract] OR "Fluorescent protein-based"[Title/Abstract] OR "Bioluminescent"[Title/Abstract] OR "Chemogenetic"[Title/Abstract])'
    function_terms = '("Sensor"[Title/Abstract] OR "Indicator"[Title/Abstract] OR "Probe"[Title/Abstract] OR "Reporter"[Title/Abstract] OR "Biosensor"[Title/Abstract])'
    specific_families = '("GCaMP" OR "GECI" OR "GEVI" OR "iSnFR" OR "GRAB" OR "dLight" OR "FRET biosensor" OR "BRET biosensor" OR "cpGFP" OR "CaBLAM")'
    noise_filter = 'NOT ("Wastewater" OR "Pollutant" OR "Review"[Publication Type])'

    # 使用 edat (录入日) 抓取最新
    search_query = f'(({core_tech} AND {function_terms}) OR {specific_families}) {noise_filter}'

    try:
        handle = Entrez.esearch(db="pubmed", term=search_query, reldate=days_back, datetype="edat", retmax=200)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list: return []

        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        papers_data = Entrez.read(handle)
        handle.close()
        
        results = []
        for paper in papers_data['PubmedArticle']:
            article = paper['MedlineCitation']['Article']
            
            # 提取基础信息
            title = article.get('ArticleTitle', 'No Title')
            journal = article.get('Journal', {}).get('Title', 'Unknown')
            pub_date = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', '2025')
            
            abstract_list = article.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            if not abstract: abstract = "[No Abstract Available Yet]"

            doi = "No DOI"
            for eid in article.get('ELocationID', []):
                if eid.attributes.get('EIdType') == 'doi':
                    doi = str(eid)

            # 统一字段名 (Lower case keys for consistency)
            results.append({
                "title": title,
                "journal": journal,
                "date": pub_date,
                "doi": f"https://doi.org/{doi}",
                "abstract": abstract
            })
            
        return results

    except Exception as e:
        print(f"⚠️ Fetch Error: {e}")
        return []
    


# ... (上面是你原有的 fetch_broad_probe_papers) ...

def fetch_papers_by_date_range(start_date, end_date):
    """
    手动模式：抓取指定日期范围的论文
    格式: start_date="2023/01/01", end_date="2023/12/31"
    """
    # === 1. 复用核心检索词 (保持与 CI 逻辑完全一致) ===
    core_tech = '("Genetically encoded"[Title/Abstract] OR "Fluorescent protein"[Title/Abstract] OR "Fluorescent protein-based"[Title/Abstract] OR "Bioluminescent"[Title/Abstract] OR "Chemogenetic"[Title/Abstract])'
    function_terms = '("Sensor"[Title/Abstract] OR "Indicator"[Title/Abstract] OR "Probe"[Title/Abstract] OR "Reporter"[Title/Abstract] OR "Biosensor"[Title/Abstract])'
    specific_families = '("GCaMP" OR "GECI" OR "GEVI" OR "iSnFR" OR "GRAB" OR "dLight" OR "FRET biosensor" OR "BRET biosensor" OR "cpGFP" OR "CaBLAM")'
    noise_filter = 'NOT ("Wastewater" OR "Pollutant" OR "Review"[Publication Type])'

    # === 2. 构建日期搜索 (Date - Entrez: 录入日期) ===
    # 相比 CI 的 reldate，这里使用具体的时间段
    date_query = f'("{start_date}"[Date - Entrez] : "{end_date}"[Date - Entrez])'
    
    search_query = f'(({core_tech} AND {function_terms}) OR {specific_families}) {noise_filter} AND {date_query}'
    
    print(f"🔍 Searching Query: {search_query}...")

    try:
        # retmax 设大一点，因为一年可能有很多篇
        handle = Entrez.esearch(db="pubmed", term=search_query, retmax=2000, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        print(f"📦 Found {len(id_list)} IDs in range {start_date} - {end_date}")
        
        if not id_list: return []

        # 分批次获取详情 (防止 URL 过长报错)
        batch_size = 200
        results = []
        
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i+batch_size]
            print(f"   Downloading batch {i}-{i+len(batch_ids)}...")
            
            handle = Entrez.efetch(db="pubmed", id=batch_ids, rettype="medline", retmode="xml")
            papers_data = Entrez.read(handle)
            handle.close()
            
            for paper in papers_data['PubmedArticle']:
                try:
                    article = paper['MedlineCitation']['Article']
                    
                    title = article.get('ArticleTitle', 'No Title')
                    journal = article.get('Journal', {}).get('Title', 'Unknown')
                    pub_date = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', 'N/A')
                    
                    abstract_list = article.get('Abstract', {}).get('AbstractText', [])
                    abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
                    if not abstract: abstract = "[No Abstract]"

                    doi = "No DOI"
                    for eid in article.get('ELocationID', []):
                        if eid.attributes.get('EIdType') == 'doi':
                            doi = str(eid)

                    results.append({
                        "title": title,
                        "journal": journal,
                        "date": pub_date,
                        "doi": f"https://doi.org/{doi}",
                        "abstract": abstract
                    })
                except Exception as e:
                    continue
                    
        return results

    except Exception as e:
        print(f"⚠️ Fetch Error: {e}")
        return []