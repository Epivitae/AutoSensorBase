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

def fetch_metadata(doi):
    """带有多重冗余的元数据抓取函数"""
    final_doi = doi.strip()
    # 提取纯数字 DOI，防止用户输入完整的 URL 导致搜索失败
    clean_doi = final_doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
    if not final_doi.startswith("http"): final_doi = f"https://doi.org/{clean_doi}"

    # --- 尝试 1: PubMed ---
    try:
        search_handle = Entrez.esearch(db="pubmed", term=clean_doi, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_record.get("IdList", [])
        if id_list:
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
        print(f"⚠️ PubMed 查询异常或未收录: {e}")

    # --- 尝试 2: Crossref API (无缝切换) ---
    print("🔄 PubMed未收录，自动切换到 Crossref API 抓取...")
    try:
        url = f"https://api.crossref.org/works/{clean_doi}"
        # 附带邮箱符合 Crossref 的 polite pool 规范，抓取速度更快且不易被封
        headers = {"User-Agent": "AutoSensorBase/1.0 (mailto:wangk@ion.ac.cn)"}
        r = requests.get(url, headers=headers, timeout=10)
        
        if r.status_code == 200:
            item = r.json().get("message", {})
            title = item.get("title", ["No Title"])[0]
            journal = item.get("container-title", ["Unknown Journal"])[0]
            
            # Crossref 时间字段较复杂，提取年份
            published = item.get("published-print", item.get("published-online", {}))
            date_parts = published.get("date-parts", [[str(datetime.datetime.now().year)]])
            year = str(date_parts[0][0])
            
            abstract = item.get("abstract", "No abstract available.")
            abstract = re.sub(r'<[^>]+>', '', abstract) # 粗略清理可能带有的 JATS XML 标签
            
            return {
                "title": title,
                "journal": journal,
                "date": year,
                "abstract": abstract,
                "doi": final_doi
            }, None
    except Exception as e:
        print(f"⚠️ Crossref 查询异常: {e}")

    # 都失败了返回错误
    return None, "所有数据库均未找到该文献"

def ai_extract(abstract):
    if not ZHIPU_API_KEY or abstract == "No abstract available.": return {}
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
        if val not in ["_No response_", "由 AI 自动提取", "None", ""]:
            return val
    return None

def main():
    body = os.environ.get("ISSUE_BODY", "")
    
    # 1. 解析基础信息
    probe_name = extract_field(body, "探针名称")
    doi = extract_field(body, "论文 DOI")
    
    if not probe_name or not doi:
        print("❌ 解析 Issue 失败：必须提供探针名称和 DOI")
        return
        
    print(f"⏳ 正在处理: {probe_name} ({doi})")
    
    # 2. 提取人工选填参数
    manual_target = extract_field(body, "检测底物")
    manual_type = extract_field(body, "探针类型")
    manual_color = extract_field(body, "发光颜色")
    manual_color_custom = extract_field(body, "其他颜色")
    
    if manual_color and "其他" in manual_color:
        manual_color = manual_color_custom
    
    # 3. 冗余抓取文献元数据
    meta, err = fetch_metadata(doi)
    
    # --- 尝试 3: 终极降级方案 ---
    if err or not meta:
        print(f"🚨 {err}。启动极客降级方案，强制创建基础条目。")
        meta = {
            "title": "Manual Verification Required",
            "journal": "Unknown",
            "date": str(datetime.datetime.now().year),
            "abstract": "No abstract available.",
            "doi": doi if doi.startswith("http") else f"https://doi.org/{doi.strip()}"
        }
        
    # 4. AI 提取
    print("🧠 正在请求 AI 提取...")
    ai_attrs = ai_extract(meta['abstract'])
    
    # 5. 智能合并
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
    
    # 6. 写入 JSON
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