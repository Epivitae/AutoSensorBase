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
    1. NO FABRICATION: If the specific probe name, target, or color is NOT explicitly written in the text, you MUST output null for that field. Do not guess.
    2. REJECT APPLICATIONS: If the authors are merely USING an existing probe for physiological/biological experiments (e.g., "we utilized GCaMP", "we expressed..."), set "is_new" to false.
    3. REJECT NON-GEFS: Strictly reject chemical dyes, PET/CT radiotracers (e.g., PSMA), nanoparticle sensors, or physical electrodes. Set "is_new" to false.

    ### INPUT
    Title: {title}
    Abstract: {abstract}

    ### OUTPUT FORMAT
    Provide a SINGLE JSON object exactly matching this structure. 
    {
    "is_new": boolean, // true ONLY if they developed/engineered a new GEFS
    "reasoning": string, // Max 15 words. Quote a short phrase proving development OR explaining why it's rejected.
    "probe_name": string | null, // Exact name (e.g., "PinkyCaMP"). If multiple, use the primary one.
    "target": string | null, // e.g., "Calcium", "Dopamine". 
    "color": string | null, // e.g., "Red", "Green". Null if unspecified.
    "type": string | null // e.g., "Ion Sensor", "Metabolite Sensor", "Voltage Sensor".
    }

    ### EXAMPLES (Study these classifications carefully)

    [Example 1: True Positive - Development (Accept)]
    Title: "A high-performance green fluorescent sensor for glutamate."
    Abstract: "We engineered iGluSnFR3, a new variant with enhanced kinetics via targeted mutagenesis..."
    Output: {
    "is_new": true,
    "reasoning": "Engineered a new variant iGluSnFR3 via targeted mutagenesis.",
    "probe_name": "iGluSnFR3",
    "target": "Glutamate",
    "color": "Green",
    "type": "Metabolite Sensor"
    }

    [Example 2: Hard Negative - Application (Must Reject)]
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

    [Example 3: Irrelevant / Out of Scope (Must Reject)]
    Title: "Synthesis of novel rhodamine derivatives for cell staining."
    Abstract: "We synthesized a new chemical dye capable of crossing the plasma membrane..."
    Output: {
    "is_new": false,
    "reasoning": "Synthesized a chemical dye, not a genetically encoded sensor.",
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