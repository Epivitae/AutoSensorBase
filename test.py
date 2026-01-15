import json
import os
import time
import toml  # 引入 toml 读取库
from data_analyzer import analyze_one_paper

RAW_FILE = "raw_papers.json"
TEST_LIMIT = 50 

# ================= 自动加载 Key (新增) =================
if not os.environ.get("ZHIPU_API_KEY"):
    try:
        secret_path = ".streamlit/secrets.toml"
        if os.path.exists(secret_path):
            with open(secret_path, "r", encoding="utf-8") as f:
                config = toml.load(f)
                if "ZHIPU_API_KEY" in config:
                    os.environ["ZHIPU_API_KEY"] = config["ZHIPU_API_KEY"]
                    print("✅ 已从 secrets.toml 加载 API Key")
    except Exception as e:
        print(f"⚠️ 无法加载 secrets.toml: {e}")
# ====================================================

def main():
    print(f"🔬 准备进行 {TEST_LIMIT} 篇的小规模测试...\n")

    if not os.path.exists(RAW_FILE):
        print("❌ 没找到 raw_papers.json")
        return

    with open(RAW_FILE, "r", encoding="utf-8") as f:
        all_papers = json.load(f)

    # 挑选前 50 篇
    target_batch = all_papers[:TEST_LIMIT]
    
    print(f"📦 选中了 {len(target_batch)} 篇文献进行测试。")
    print("=" * 60)

    new_probe_count = 0
    
    for i, paper in enumerate(target_batch):
        title = paper.get('title', 'No Title')[:60]
        print(f"[{i+1}/{TEST_LIMIT}] 分析中: {title}...")

        try:
            result = analyze_one_paper(paper)

            if result:
                is_new = result.get('is_new')
                reasoning = result.get('reasoning', 'No reasoning provided')
                probe_name = result.get('probe_name', 'N/A')

                if is_new:
                    print(f"   🎉 [判定: NEW] | 探针: {probe_name}")
                    print(f"   🧠 推理: \033[92m{reasoning}\033[0m") 
                    new_probe_count += 1
                else:
                    print(f"   ❌ [判定: REJECT]")
                    print(f"   🧠 推理: \033[90m{reasoning}\033[0m") 
            else:
                print("   ⚠️ AI 响应为空 (可能是 Key 依然有问题)")

        except Exception as e:
            print(f"   ⚠️ 出错: {e}")
        
        print("-" * 60)
        time.sleep(0.5)

    print(f"\n📊 测试结束！")
    print(f"   样本数: {len(target_batch)}")
    print(f"   新探针数: {new_probe_count}")

if __name__ == "__main__":
    main()