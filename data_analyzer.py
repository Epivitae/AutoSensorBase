import json
import requests
import os

# 锁定稳定版本模型，保障 GitHub Actions 自动化长治久安
MODEL_NAME = "gemini-2.0-flash-001" 

def analyze_one_paper(paper):
    # ✅ 每次调用函数时实时获取 Gemini Key
    api_key = os.environ.get("GEMINI_API_KEY")

    if not api_key:
        print("❌ Error: GEMINI_API_KEY not found in environment variables.")
        return None

    # 🌟 构造 Gemini 官方 REST URL
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{MODEL_NAME}:generateContent?key={api_key}"
    
    headers = {
        "Content-Type": "application/json"
    }
    
    # 🌟 终极版 Prompt (保持你的原版优质 Prompt 不变)
    prompt = f"""
    You are a strict Senior Editor at *Nature Methods*. Your job is to filter papers to find ONLY those describing the **Development of NEW Genetically Encoded Fluorescent Sensors**.
    
    --- INPUT DATA ---
    Title: {paper.get('title', '')}
    Abstract: {paper.get('abstract', '')}
    
    --- YOUR REASONING STEPS (Mental Sandbox) ---
    1. **Check Material**: Is it a protein? (Exclude: Chemical dyes, Nanoparticles, Electrodes).
    2. **Check Novelty**: Is it a NEW variant created by the authors? (Exclude: Existing sensors like GCaMP6, dLight1 used for experiments).
    3. **Check Action**: Did they perform "Engineering", "Screening", "Mutagenesis", or "Directed Evolution"? (Exclude: "Used", "Applied", "Validated").
    4. **Check Name**: Did they name a new sensor? (e.g., "GRAB-DA3.0"). If they just say "A sensor", be suspicious.
    
    --- OUTPUT REQUIREMENTS ---
    Return a SINGLE JSON object. No Markdown. No text outside JSON.
    
    The JSON must contain:
    - "reasoning": A short sentence explaining your judgment.
    - "evidence_quote": A short quote from the text proving development (or "N/A").
    - "is_new": boolean (true/false).
    - "probe_name": String (The specific NEW name, or "Unknown").
    - "target": String.
    - "color": String.
    - "type": String (e.g., "Intensometric", "FRET", "Bioluminescent").

    --- EXAMPLES ---
    [Case 1: Application - REJECT]
    Title: "Imaging dopamine in striatum using dLight1."
    Abstract: "We used the dLight1 sensor to monitor..."
    Result: {{"reasoning": "Authors used an existing sensor dLight1, no new engineering.", "is_new": false}}

    [Case 2: Development - ACCEPT]
    Title: "A sensitive red fluorescent calcium indicator."
    Abstract: "We performed site-directed mutagenesis and screened 1000 variants to create jRGECO1a..."
    Result: {{"reasoning": "Authors performed mutagenesis and screening to create a named variant jRGECO1a.", "is_new": true, "probe_name": "jRGECO1a", "target": "Calcium", "color": "Red", "type": "Intensometric"}}

    [Case 3: No Abstract - TITLE JUDGMENT]
    Title: "Genetically encoded sensors for serotonin."
    Abstract: "[No Abstract]"
    Result: {{"reasoning": "Title strongly implies a new class of sensors.", "is_new": true, "probe_name": "Serotonin Sensor (Unnamed)", "target": "Serotonin", "color": "Unknown", "type": "Genetically Encoded Sensor"}}
    """
    
    # 🌟 适配 Gemini 的 Payload 结构
    payload = {
        "contents": [{
            "parts": [{"text": prompt}]
        }],
        "generationConfig": {
            "responseMimeType": "application/json", # 强制只输出纯 JSON，拒绝 Markdown 包装
            "temperature": 0.01,
            "topP": 0.1
        }
    }

    try:
        resp = requests.post(url, headers=headers, json=payload, timeout=30)
        
        if resp.status_code == 200:
            # 🌟 解析 Gemini 的返回体路径
            content = resp.json()['candidates'][0]['content']['parts'][0]['text']
            
            try:
                # 不再需要 replace("```json")，Gemini 在配置后会直接返回干净的 JSON 字符串
                result = json.loads(content.strip())
                
                # --- Guardrail (保持你的后处理逻辑) ---
                if result.get('is_new'):
                    reasoning = result.get('reasoning', '').lower()
                    if "used" in reasoning and "existing" in reasoning and "develop" not in reasoning:
                        print(f"   🛡️ [Guardrail] AI hallucinated 'True' but reasoning implies 'False'. Correcting...")
                        result['is_new'] = False
                
                return result
            except json.JSONDecodeError:
                print(f"⚠️ JSON Parse Error: {content[:50]}...")
                return None
        else:
            print(f"⚠️ API Error: HTTP {resp.status_code} - {resp.text}")
            
    except Exception as e:
        # 打印错误时，用 replace 把真实 Key 隐藏掉
        safe_url = url.replace(api_key, "HIDDEN_KEY") if api_key else url
        print(f"⚠️ Network Error (URL: {safe_url}): {e}")
    
    return None