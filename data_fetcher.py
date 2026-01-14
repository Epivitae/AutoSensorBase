import time
from Bio import Entrez
import datetime

# 设置邮箱（PubMed 要求）
Entrez.email = "wangk@ion.ac.cn"

def fetch_recent_papers(days_back=7):
    """
    爬取过去 N 天的论文
    """
    # 动态生成日期范围
    today = datetime.date.today()
    start_date = today - datetime.timedelta(days=days_back)
    date_range = f'("{start_date.strftime("%Y/%m/%d")}"[Date - Publication] : "{today.strftime("%Y/%m/%d")}"[Date - Publication])'

    # 搜索词 (优化版)
    search_term = (
        '("genetically encoded"[Title/Abstract] AND '
        '("sensor"[Title/Abstract] OR "indicator"[Title/Abstract] OR "probe"[Title/Abstract])) '
        'AND ("develop*"[Title/Abstract] OR "engineer*"[Title/Abstract] OR "design*"[Title/Abstract] OR "new"[Title/Abstract]) '
        f'AND {date_range}'
    )
    
    print(f"🕷️ Search Term: {search_term}")
    
    try:
        # 1. 搜索
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100, sort='date')
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        
        if not id_list:
            return []

        # 2. 获取详情
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        papers_data = Entrez.read(handle)
        handle.close()
        
        results = []
        for paper in papers_data['PubmedArticle']:
            article = paper['MedlineCitation']['Article']
            
            # 提取 DOI
            doi = "No DOI"
            for eid in article.get('ELocationID', []):
                if eid.attributes.get('EIdType') == 'doi':
                    doi = str(eid)
            
            # 提取摘要
            abstract_list = article.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)

            if len(abstract) < 50: continue # 跳过无摘要的

            results.append({
                "title": article.get('ArticleTitle', 'No Title'),
                "abstract": abstract,
                "doi": f"https://doi.org/{doi}",
                "journal": article.get('Journal', {}).get('Title', 'Unknown'),
                "date": article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', '2025')
            })
            
        return results

    except Exception as e:
        print(f"⚠️ Error fetching papers: {e}")
        return []