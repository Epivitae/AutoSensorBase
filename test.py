import os
import json
import requests

# 1. 强行给 Python 挂上你配置好的本地代理隧道（保持你的代理软件开启）
os.environ["HTTP_PROXY"] = "http://127.0.0.1:7890"  # 如果你的端口不是 7890，请改成你的真实端口
os.environ["HTTPS_PROXY"] = "http://127.0.0.1:7890"

# 2. 换上你刚申请的、纯净的 AI Studio Gemini Key
os.environ["GEMINI_API_KEY"] = "AIzaSyD-utT404JxrO4Me6LP35fp6Qrx7MEoj9Y"

api_key = os.environ.get("GEMINI_API_KEY")
print(f"📡 正在检测本地 Gemini API Key...")
if not api_key or "你的新Key" in api_key:
    print("❌ 错误：请先在代码中填入真实的 GEMINI_API_KEY！")
    exit(1)

# 模拟摘要
mock_abstract = (
    "Genetically encoded calcium indicators (GECIs) are critical tools for imaging neural activity. "
    "Here we develop GCaMP9, a novel green fluorescent protein-based sensor with ultra-high sensitivity "
    "and faster kinetics for monitoring intracellular Calcium dynamics in vivo."
)

def test_gemini_fixed(abstract):
    # 🌟 严格对齐官方的标准 REST URL 格式
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={api_key}"
    
    headers = {
        "Content-Type": "application/json"
    }
    
    prompt = f"""Read this abstract and extract the fluorescent sensor's attributes.
    Abstract: {abstract}
    
    Return JSON ONLY with these keys:
    {{
        "target": "What molecule is detected? (e.g. Dopamine, Calcium)",
        "color": "Fluorescence color (Green, Red, Blue, Yellow, etc.)",
        "type": "Sensor type (e.g. CPFP, GPCR-based, snifit)"
    }}
    If unsure, use "Unknown". Do NOT return Markdown block (do not include ```json)."""
    
    payload = {
        "contents": [{
            "parts": [{"text": prompt}]
        }],
        "generationConfig": {
            "responseMimeType": "application/json", # 强行约束返回纯 JSON
            "temperature": 0.1
        }
    }

    print("\n🧠 正在向 Google Gemini 官方标准端点发送请求...")
    try:
        resp = requests.post(url, headers=headers, json=payload, timeout=15)
        print(f"📬 服务器响应状态码: {resp.status_code}")
        
        if resp.status_code == 200:
            res_json = resp.json()
            text_content = res_json['candidates'][0]['content']['parts'][0]['text']
            return json.loads(text_content.strip()), None
        else:
            return None, f"HTTP {resp.status_code}: {resp.text}"
    except Exception as e:
        return None, f"连接异常: {str(e)}"

# 执行
result, error = test_gemini_fixed(mock_abstract)
print("\n--- 📝 Gemini 测试报告 ---")
if error:
    print(f"❌ 测试失败！原因：{error}")
else:
    print("🎉 恭喜！Gemini 官方接口完全激活跑通！解析结果如下：")
    print(json.dumps(result, indent=4, ensure_ascii=False))