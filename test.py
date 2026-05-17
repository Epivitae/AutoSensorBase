import json
import requests
import os
import time

# 1. 强行挂上本地代理（排查嫌疑人一：本地直连被拦截）
os.environ["HTTP_PROXY"] = "http://127.0.0.1:7890"  # 确保这里是你的代理软件真实端口
os.environ["HTTPS_PROXY"] = "http://127.0.0.1:7890"

# 2. 填入你刚刚新申请的 Key（排查嫌疑人二：旧 Key 被关小黑屋）
os.environ["GEMINI_API_KEY"] = "AIzaSyAHOlPMRLwZb0uVV1iVhy73xZf6nSyWUWM"

# 3. 使用最基础的通用别名（排查嫌疑人三：带 001 后缀的模型额度 Bug）
MODEL_NAME = "gemini-2.0-flash" 

def test_gemini_connection():
    api_key = os.environ.get("GEMINI_API_KEY")
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{MODEL_NAME}:generateContent?key={api_key}"
    headers = {"Content-Type": "application/json"}
    
    # 构造一条极简的测试数据
    mock_paper = {
        "title": "A sensitive red fluorescent calcium indicator.",
        "abstract": "We performed site-directed mutagenesis and screened 1000 variants to create jRGECO1a, a new sensor for in vivo imaging."
    }
    
    prompt = f"""
    You are a strict Senior Editor at *Nature Methods*. Your job is to filter papers to find ONLY those describing the **Development of NEW Genetically Encoded Fluorescent Sensors**.
    
    Title: {mock_paper['title']}
    Abstract: {mock_paper['abstract']}
    
    Return a SINGLE JSON object with keys: "reasoning", "is_new" (boolean), "probe_name", "target", "color", "type".
    """
    
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "responseMimeType": "application/json",
            "temperature": 0.01
        }
    }

    print("🚀 开始诊断：正在连接 Google Gemini 服务器...")
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            resp = requests.post(url, headers=headers, json=payload, timeout=15)
            print(f"📬 收到服务器状态码: {resp.status_code}")
            
            if resp.status_code == 200:
                print("✅ 恭喜！连接完全正常！解析结果如下：")
                content = resp.json()['candidates'][0]['content']['parts'][0]['text']
                print(json.dumps(json.loads(content.strip()), indent=4, ensure_ascii=False))
                return
                
            elif resp.status_code == 429:
                # 打印详细错误信息，帮我们诊断
                error_msg = resp.json().get('error', {}).get('message', '')
                print(f"⚠️ 触发限制 (429)。详细信息: {error_msg}")
                
                if "limit: 0" in error_msg:
                    print("🚨 致命错误：依然是 Limit 0！")
                    print("   👉 诊断结果：如果代理没问题，说明你的这个 API Key 依然处于被锁死的状态（或者该账号没开通此模型权限），请换个新号/新Key！")
                    return
                else:
                    wait_time = 5 * (attempt + 1)
                    print(f"   ⏳ 只是请求过快，等待 {wait_time} 秒后重试...")
                    time.sleep(wait_time)
            else:
                print(f"❌ 发生其他错误：{resp.text}")
                return
                
        except Exception as e:
            print(f"❌ 网络异常 (代理可能没配置对): {e}")
            return
            
    print("❌ 三次尝试全部失败。")

if __name__ == "__main__":
    test_gemini_connection()