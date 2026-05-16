import json
import requests
import os
import time  # ✅ 新增引入 time 模块，用于控制请求频率

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
    
    # 🌟 终极版 Prompt
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
            "responseMimeType": "application/json", # 强制只输出纯 JSON
            "temperature": 0.01,
            "topP": 0.1
        }
    }

    # 🌟 新增：重试与限流保护机制 (最多尝试3次)
    max_retries = 3
    for attempt in range(max_retries):
        try:
            resp = requests.post(url, headers=headers, json=payload, timeout=30)
            
            if resp.status_code == 200:
                # 解析返回结果
                content = resp.json()['candidates'][0]['content']['parts'][0]['text']
                
                try:
                    result = json.loads(content.strip())
                    
                    # --- Guardrail ---
                    if result.get('is_new'):
                        reasoning = result.get('reasoning', '').lower()
                        if "used" in reasoning and "existing" in reasoning and "develop" not in reasoning:
                            print(f"   🛡️ [Guardrail] AI hallucinated 'True' but reasoning implies 'False'. Correcting...")
                            result['is_new'] = False
                    
                    # ✅ 成功获取结果，强制休眠 4.5 秒（保证一分钟最多 13~14 次，避开 15 次的限制）
                    time.sleep(4.5)
                    return result
                    
                except json.JSONDecodeError:
                    print(f"⚠️ JSON Parse Error: {content[:50]}...")
                    time.sleep(4.5)
                    return None
                    
            elif resp.status_code == 429:
                # 触发了频率限制
                wait_time = 5 * (attempt + 1) # 渐进式等待：5秒、10秒、15秒
                print(f"   ⏳ 触发频率限制 (HTTP 429)，等待 {wait_time} 秒后重试 ({attempt + 1}/{max_retries})...")
                time.sleep(wait_time)
                continue # 继续下一次循环重试
                
            else:
                # 其他 API 报错 (如 400, 500 等)
                print(f"⚠️ API Error: HTTP {resp.status_code} - {resp.text}")
                time.sleep(4.5)
                return None
                
        except Exception as e:
            # 网络断开或超时报错
            safe_url = url.replace(api_key, "HIDDEN_KEY") if api_key else url
            print(f"⚠️ Network Error (URL: {safe_url}): {e}")
            wait_time = 5 * (attempt + 1)
            print(f"   ⏳ 网络异常，等待 {wait_time} 秒后重试...")
            time.sleep(wait_time)

    # 如果三次都失败了，跳出循环，返回 None
    print("   ❌ 多次请求失败，放弃该篇文献。")
    time.sleep(4.5)
    return None