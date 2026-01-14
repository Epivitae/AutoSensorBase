import json
import requests
import os

# 从环境变量读取 Key，确保安全
API_KEY = os.environ.get("ZHIPU_API_KEY") 
API_URL = "https://open.bigmodel.cn/api/paas/v4/chat/completions"
MODEL_NAME = "glm-4-flash"

def analyze_one_paper(paper):
    if not API_KEY:
        print("❌ Error: ZHIPU_API_KEY not found in environment variables.")
        return None

    headers = {"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"}
    
    prompt = f"""
    Title: {paper['title']}
    Abstract: {paper['abstract']}
    
    Task: Did they DEVELOP a NEW genetically encoded fluorescent sensor?
    (Ignore reviews and pure applications).
    
    If YES, return JSON:
    {{
        "is_new": true,
        "probe_name": "Name",
        "target": "Molecule",
        "color": "Color",
        "type": "Type"
    }}
    
    If NO, return: {{"is_new": false}}
    
    Return ONLY JSON string. No Markdown.
    """
    
    payload = {
        "model": MODEL_NAME,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.1
    }

    try:
        resp = requests.post(API_URL, headers=headers, json=payload, timeout=30)
        if resp.status_code == 200:
            content = resp.json()['choices'][0]['message']['content']
            clean_json = content.replace("```json", "").replace("```", "").strip()
            return json.loads(clean_json)
    except Exception as e:
        print(f"⚠️ AI Analysis Failed: {e}")
    
    return None