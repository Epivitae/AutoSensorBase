import json
import requests
import os

# 不要在这里读 Key，否则 import 的时候就定死了
# API_KEY = os.environ.get("ZHIPU_API_KEY")  <-- 删掉这行

API_URL = "https://open.bigmodel.cn/api/paas/v4/chat/completions"
MODEL_NAME = "glm-4-flash" 

def analyze_one_paper(paper):
    # ✅ 改到这里读取！每次调用函数时实时获取 Key
    api_key = os.environ.get("ZHIPU_API_KEY")

    if not api_key:
        print("❌ Error: ZHIPU_API_KEY not found in environment variables.")
        return None

    headers = {
        "Authorization": f"Bearer {api_key}", 
        "Content-Type": "application/json"
    }
    
    # 🌟 终极版 Prompt (保持不变)
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
    
    payload = {
        "model": MODEL_NAME,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.01,
        "top_p": 0.1
    }

    try:
        resp = requests.post(API_URL, headers=headers, json=payload, timeout=30)
        
        if resp.status_code == 200:
            content = resp.json()['choices'][0]['message']['content']
            clean_content = content.replace("```json", "").replace("```", "").strip()
            
            try:
                result = json.loads(clean_content)
                
                # --- Guardrail ---
                if result.get('is_new'):
                    reasoning = result.get('reasoning', '').lower()
                    if "used" in reasoning and "existing" in reasoning and "develop" not in reasoning:
                        print(f"   🛡️ [Guardrail] AI hallucinated 'True' but reasoning implies 'False'. Correcting...")
                        result['is_new'] = False
                
                return result
            except json.JSONDecodeError:
                print(f"⚠️ JSON Parse Error: {clean_content[:50]}...")
                return None
        else:
            print(f"⚠️ API Error: {resp.status_code}")
            
    except Exception as e:
        print(f"⚠️ Network Error: {e}")
    
    return None