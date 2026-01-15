import time
from Bio import Entrez
import datetime

# 设置邮箱
Entrez.email = "wangk@ion.ac.cn"

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