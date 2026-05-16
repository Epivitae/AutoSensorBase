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
    if not ZHIPU_API_KEY: return {}
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
    return {}

def extract_field(body, keyword):
    """从 Issue 正文中提取特定字段的值"""
    match = re.search(rf'### {keyword}[^\n]*\n+([^\n]+)', body)
    if match:
        val = match.group(1).strip()
        # 过滤掉未填写的默认值
        if val not in ["_No response_", "由 AI 自动提取", "None", ""]:
            return val
    return None

def main():
    body = os.environ.get("ISSUE_BODY", "")
    
    # 解析 Issue 表单内容
    probe_name = extract_field(body, "探针名称")
    doi = extract_field(body, "论文 DOI")
    
    if not probe_name or not doi:
        print("❌ 解析 Issue 失败：必须提供探针名称和 DOI")
        return
        
    print(f"⏳ 正在处理: {probe_name} ({doi})")
    
    # 获取人工选填的参数
    manual_target = extract_field(body, "检测底物")
    manual_color = extract_field(body, "发光颜色")
    manual_type = extract_field(body, "探针类型")
    
    # 抓取文献元数据
    meta, err = fetch_pubmed(doi)
    if err:
        print(f"❌ 抓取失败: {err}")
        return
        
    # AI 提取
    print("🧠 正在请求 AI 提取...")
    ai_attrs = ai_extract(meta['abstract'])
    
    # 智能合并（人工填写 > AI 提取 > 默认值 Unknown）
    final_target = manual_target if manual_target else ai_attrs.get("target", "Unknown")
    final_color = manual_color if manual_color else ai_attrs.get("color", "Unknown")
    final_type = manual_type if manual_type else ai_attrs.get("type", "Unknown")
    
    new_entry = {
        "probe_name": probe_name,
        "is_new": True,
        **meta,
        "target": final_target,
        "color": final_color,
        "type": final_type
    }
    
    # 更新 JSON
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r", encoding="utf-8") as f: data = json.load(f)
    else:
        data = []
        
    data.insert(0, new_entry)
    
    with open(DATA_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
    print("✅ 探针数据合并与更新成功！")

if __name__ == "__main__":
    main()