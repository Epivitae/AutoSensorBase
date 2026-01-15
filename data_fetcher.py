import time
from Bio import Entrez
import datetime

# 设置邮箱
Entrez.email = "wangk@ion.ac.cn"

# ================= 1. 旧版函数 (保持不变用于对比) =================
def fetch_recent_papers(days_back=7):
    # ... (保持你原来旧版代码的内容，为了debug_compare能运行，这里就不重复贴了，保留你文件里的原样即可) ...
    # 如果你刚才覆盖了，请确保这里有旧代码，或者直接看下面的新函数替换掉原来的新函数
    
    # 为了方便，我这里只贴出你需要替换的【新版函数】
    pass 

# ================= 2. 新版函数 (V3 修复版) =================
def fetch_broad_probe_papers(days_back=7):
    """
    V3 改进版：修复生物发光漏抓，增加环境噪声过滤
    """
    
    # 1. 核心技术 (扩充：涵盖生物发光、化学遗传、光遗传)
    core_tech = '(' \
                '"Genetically encoded"[Title/Abstract] OR ' \
                '"Fluorescent protein"[Title/Abstract] OR ' \
                '"Fluorescent protein-based"[Title/Abstract] OR ' \
                '"Bioluminescent"[Title/Abstract] OR ' \
                '"Bioluminescence"[Title/Abstract] OR ' \
                '"Chemogenetic"[Title/Abstract] OR ' \
                '"Optogenetic"[Title/Abstract]' \
                ')'
    
    # 2. 核心功能 (扩充)
    function_terms = '(' \
                     '"Sensor"[Title/Abstract] OR ' \
                     '"Indicator"[Title/Abstract] OR ' \
                     '"Probe"[Title/Abstract] OR ' \
                     '"Reporter"[Title/Abstract] OR ' \
                     '"Biosensor"[Title/Abstract] OR ' \
                     '"Integrator"[Title/Abstract]' \
                     ')'
    
    # 3. 具体家族 (黑话列表，持续更新)
    specific_families = '(' \
                        '"GCaMP"[Title/Abstract] OR ' \
                        '"GECI"[Title/Abstract] OR ' \
                        '"GEVI"[Title/Abstract] OR ' \
                        '"iSnFR"[Title/Abstract] OR ' \
                        '"GRAB"[Title/Abstract] OR ' \
                        '"dLight"[Title/Abstract] OR ' \
                        '"FRET biosensor"[Title/Abstract] OR ' \
                        '"BRET biosensor"[Title/Abstract] OR ' \
                        '"cpGFP"[Title/Abstract] OR ' \
                        '"CaBLAM"[Title/Abstract]' \
                        ')'

    # 4. 噪声过滤 (新增：剔除环境监测、纳米材料等非遗传编码领域)
    # 注意：不要过滤 "Virus"，因为很多探针用 Virus 载体
    noise_filter = 'NOT (' \
                   '"Wastewater"[Title/Abstract] OR ' \
                   '"Pollutant"[Title/Abstract] OR ' \
                   '"Electrochemical"[Title/Abstract] OR ' \
                   '"Nanoparticle"[Title/Abstract] OR ' \
                   '"Polymer"[Title/Abstract] OR ' \
                   '"Review"[Publication Type]' \
                   ')'

    # 组合检索词
    search_query = f'(({core_tech} AND {function_terms}) OR {specific_families}) {noise_filter}'
    
    # 本地打标关键词
    dev_keywords = ["develop", "engineer", "design", "novel", "new", "variant", "optimize", "screen", "characteriz"]

    try:
        # 使用 edat (录入日)
        handle = Entrez.esearch(db="pubmed", term=search_query, reldate=days_back, datetype="edat", retmax=500)
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
            
            title = article.get('ArticleTitle', 'No Title')
            journal = article.get('Journal', {}).get('Title', 'Unknown')
            
            abstract_list = article.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            
            doi = "No DOI"
            for eid in article.get('ELocationID', []):
                if eid.attributes.get('EIdType') == 'doi':
                    doi = str(eid)

            # 智能打标
            text_for_search = (title + " " + abstract).lower()
            is_dev_paper = any(kw in text_for_search for kw in dev_keywords)
            
            if not abstract: abstract = "[No Abstract Available Yet]"

            results.append({
                "Title": title,
                "title": title, 
                "Journal": journal,
                "DOI": f"https://doi.org/{doi}",
                "doi": f"https://doi.org/{doi}", 
                "Is_Development_Likely": "✅" if is_dev_paper else "",
                "Abstract": abstract,
                "abstract": abstract 
            })
            
        return results

    except Exception as e:
        print(f"⚠️ New Fetch Error: {e}")
        return []

# 请确保把 fetch_recent_papers 也保留在文件里，不要删，否则 debug_compare.py 会报错