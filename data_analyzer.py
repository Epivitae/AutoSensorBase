import json
import requests
import os
import time

# 锁定 DeepSeek 的对话模型
MODEL_NAME = "deepseek-chat" 

def analyze_one_paper(paper):
    # ✅ 获取 DeepSeek Key
    api_key = os.environ.get("DEEPSEEK_API_KEY")

    if not api_key:
        print("❌ Error: DEEPSEEK_API_KEY not found in environment variables.")
        return None

    # 🌟 修改点 1：DeepSeek 的官方 API 端点
    url = "https://api.deepseek.com/chat/completions"
    
    # 🌟 修改点 2：DeepSeek 鉴权必须在 Header 里使用 Bearer Token
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    
    # 🌟 将 Prompt 拆分为 System 和 User，这在 DeepSeek/OpenAI 架构下效果更好
    system_prompt = """
    You are an expert bioimage analyst. Your objective is to extract structured data ONLY if the paper describes the ORIGINAL DEVELOPMENT, ENGINEERING, or DIRECT OPTIMIZATION of a NEW genetically encoded fluorescent sensor (GEFS).

    ### CRITICAL RULES (Anti-Hallucination & Filtering)
    1. NO FABRICATION: If the specific probe name or target is NOT explicitly written in the text, you MUST output null for that field. Do not guess targets.
    2. REJECT APPLICATIONS: If the authors are merely USING an existing probe for physiological/biological experiments (e.g., "we utilized GCaMP", "we expressed..."), set "is_new" to false.
    3. REJECT NON-GEFS: Strictly reject chemical dyes, PET/CT radiotracers (e.g., PSMA), nanoparticle sensors, or physical electrodes. Set "is_new" to false.
    
    4. COLOR INFERENCE (CRITICAL): Abstracts often use biological jargon or specific nomenclature instead of basic color names. You MUST infer the base color category using these rules:
       - Green: Explicitly mentions GFP, EGFP, sfGFP, mNeonGreen, YFP, Venus, "green fluorescent", or emission ~500-530nm. 
         [NOMENCLATURE RULE]: If the probe name ends with suffixes like "-SnFR", "-CaMP", or starts with "i" (e.g., iGlucoSnFR, GCaMP) AND there is no indication of it being red-shifted, DEFAULT to "Green".
       - Red: Explicitly mentions RFP, mCherry, mApple, mRuby, RGECO, "red fluorescent", or emission ~580-620nm. 
         [NOMENCLATURE RULE]: Probe names starting with "R-", "jR-", or "m-" alongside a standard suffix (e.g., R-iGluSnFR, jRGECO) denote "Red".
       - Cyan: CFP, mTurquoise, emission ~470-490nm.
       - NIR: iRFP, miRFP, near-infrared, emission >650nm.
       If multiple colors are developed, list the primary one or separate by a slash (e.g., "Green/Red"). Output `null` ONLY if no keywords or nomenclature rules apply.

    ### INPUT
    Title: {title}
    Abstract: {abstract}

    ### OUTPUT FORMAT
    Provide a SINGLE JSON object exactly matching this structure. 
    {
    "is_new": boolean, // true ONLY if they developed/engineered a new GEFS
    "reasoning": string, // Max 15 words. Quote a short phrase proving development OR explaining why it's rejected.
    "probe_name": string | null, // Exact name (e.g., "PinkyCaMP"). If multiple, use the primary one.
    "target": string | null, // e.g., "Calcium", "Dopamine", "Glucose". 
    "color": string | null, // Base category: "Green", "Red", "Cyan", "Yellow", "NIR". Infer based on Rule 4. 
    "type": string | null // e.g., "Ion Sensor", "Metabolite Sensor", "Voltage Sensor", "Neurotransmitter Sensor".
    }

    ### EXAMPLES (Study these classifications carefully)

    [Example 1: Hard Case - Color Inference via Nomenclature]
    Title: "iGlucoSnFR2: A genetically encoded fluorescent sensor for measuring glucose."
    Abstract: "We have developed a second generation, genetically encoded intensity-based glucose sensing fluorescent reporter (iGlucoSnFR2)..."
    Output: {
    "is_new": true,
    "reasoning": "Developed a second generation glucose sensor iGlucoSnFR2.",
    "probe_name": "iGlucoSnFR2",
    "target": "Glucose",
    "color": "Green", 
    "type": "Metabolite Sensor"
    }

    [Example 2: True Positive - Standard Red Probe]
    Title: "A red fluorescent voltage indicator with fast kinetics."
    Abstract: "Here we report the development of JEDI-2P, a red-shifted voltage sensor with enhanced two-photon cross-section..."
    Output: {
    "is_new": true,
    "reasoning": "Report the development of JEDI-2P, a red-shifted voltage sensor.",
    "probe_name": "JEDI-2P",
    "target": "Voltage",
    "color": "Red",
    "type": "Voltage Sensor"
    }

    [Example 3: Hard Negative - Application (Must Reject)]
    Title: "In vivo imaging of neural dynamics in awake mice."
    Abstract: "We expressed the calcium indicator GCaMP6s in the visual cortex to monitor activity..."
    Output: {
    "is_new": false,
    "reasoning": "Used an existing indicator GCaMP6s for biological monitoring; no new development.",
    "probe_name": null,
    "target": null,
    "color": null,
    "type": null
    }
    """
    
    user_prompt = f"""
    --- INPUT DATA ---
    Title: {paper.get('title', '')}
    Abstract: {paper.get('abstract', '')}
    """
    
    # 🌟 修改点 3：适配 DeepSeek/OpenAI 的 Payload 结构
    payload = {
        "model": MODEL_NAME,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        "response_format": {"type": "json_object"}, # 强制只输出纯 JSON
        "temperature": 0.01,
        "top_p": 0.1
    }

    # 重试与限流保护机制 (最多尝试3次)
    max_retries = 3
    for attempt in range(max_retries):
        try:
            resp = requests.post(url, headers=headers, json=payload, timeout=40)
            
            if resp.status_code == 200:
                # 🌟 修改点 4：解析 DeepSeek 的返回结构
                content = resp.json()['choices'][0]['message']['content']
                
                try:
                    result = json.loads(content.strip())
                    
                    # --- Guardrail ---
                    if result.get('is_new'):
                        reasoning = result.get('reasoning', '').lower()
                        if "used" in reasoning and "existing" in reasoning and "develop" not in reasoning:
                            print(f"   🛡️ [Guardrail] AI hallucinated 'True' but reasoning implies 'False'. Correcting...")
                            result['is_new'] = False
                    
                    # DeepSeek 的并发限制可能和 Gemini 不同，这里保持 4.5 秒比较稳妥
                    time.sleep(4.5)
                    return result
                    
                except json.JSONDecodeError:
                    print(f"⚠️ JSON Parse Error: {content[:50]}...")
                    time.sleep(4.5)
                    return None
                    
            elif resp.status_code == 429:
                wait_time = 5 * (attempt + 1)
                print(f"   ⏳ 触发频率限制 (HTTP 429)，等待 {wait_time} 秒后重试 ({attempt + 1}/{max_retries})...")
                time.sleep(wait_time)
                continue
                
            else:
                print(f"⚠️ API Error: HTTP {resp.status_code} - {resp.text}")
                time.sleep(4.5)
                return None
                
        except Exception as e:
            # 隐藏 key 打印错误
            safe_url = url
            print(f"⚠️ Network Error: {e}")
            wait_time = 5 * (attempt + 1)
            print(f"   ⏳ 网络异常，等待 {wait_time} 秒后重试...")
            time.sleep(wait_time)

    print("   ❌ 多次请求失败，放弃该篇文献。")
    time.sleep(4.5)
    return None