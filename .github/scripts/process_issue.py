import os
import re
import json
import requests
from Bio import Entrez
import datetime

# 初始化
Entrez.email = "wangk@ion.ac.cn"
ZHIPU_API_KEY = os.environ.get("ZHIPU_API_KEY")
DATA_FILE = "processed_probes.json"

def fetch_pubmed(doi):
    try:
        final_doi = doi.strip()
        if not final_doi.startswith("http"): final_doi = f"https://doi.org/{final_doi}"
        
        search_handle = Entrez.esearch(db="pubmed", term=doi, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_record.get("IdList", [])
        if not id_list: return None, "PubMed未找到该DOI"

        fetch_handle = Entrez.efetch(db="pubmed", id=id_list[0], rettype="medline", retmode="xml")
        data = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        article = data['PubmedArticle'][0]['MedlineCitation']['Article']
        abstract_list = article.get('Abstract', {}).get('AbstractText', [])
        abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
        
        return {
            "title": article.get('ArticleTitle', 'No Title'),
            "journal": article.get('Journal', {}).get('Title', 'Unknown Journal'),
            "date": article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', str(datetime.datetime.now().year)),
            "abstract": abstract,
            "doi": final_doi
        }, None
    except Exception as e:
        return None, str(e)

def ai_extract(abstract):
    if not ZHIPU_API_KEY: return {"target": "Unknown", "color": "Unknown", "type": "Unknown"}
    try:
        url = "https://open.bigmodel.cn/api/paas/v4/chat/completions"
        headers = {"Authorization": f"Bearer {ZHIPU_API_KEY}", "Content-Type": "application/json"}
        prompt = f"""Read this abstract and extract attributes. Abstract: {abstract}
        Return JSON ONLY: {{"target": "Detected molecule", "color": "Fluorescence color", "type": "Sensor type"}}"""
        
        resp = requests.post(url, headers=headers, json={"model": "glm-4-flash", "messages": [{"role": "user", "content": prompt}], "temperature": 0.1}, timeout=15)
        if resp.status_code == 200:
            content = resp.json()['choices'][0]['message']['content']
            clean_json = content.replace("```json", "").replace("```", "").strip()
            return json.loads(clean_json)
    except:
        pass
    return {"target": "Unknown", "color": "Unknown", "type": "Unknown"}

def main():
    body = os.environ.get("ISSUE_BODY", "")
    
    # 解析 Issue 表单内容
    name_match = re.search(r'### 探针名称[^\n]*\n+([^\n]+)', body)
    doi_match = re.search(r'### 论文 DOI[^\n]*\n+([^\n]+)', body)
    
    if not name_match or not doi_match:
        print("❌ 解析 Issue 失败")
        return
        
    probe_name = name_match.group(1).strip()
    doi = doi_match.group(1).strip()
    
    print(f"⏳ 正在处理: {probe_name} ({doi})")
    
    meta, err = fetch_pubmed(doi)
    if err:
        print(f"❌ 抓取失败: {err}")
        return
        
    print("🧠 正在请求 AI 提取...")
    ai_attrs = ai_extract(meta['abstract'])
    
    new_entry = {"probe_name": probe_name, "is_new": True, **meta, **ai_attrs}
    
    # 更新 JSON
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f: data = json.load(f)
    else:
        data = []
        
    data.insert(0, new_entry)
    
    with open(DATA_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
    print("✅ JSON 更新成功！")

if __name__ == "__main__":
    main()