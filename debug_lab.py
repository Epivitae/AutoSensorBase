import os
import json
import toml
from data_analyzer import analyze_one_paper

# ================= 配置区 =================
# 自动尝试加载本地 secrets.toml，免去手动输 Key 的麻烦
if not os.environ.get("ZHIPU_API_KEY"):
    try:
        if os.path.exists(".streamlit/secrets.toml"):
            with open(".streamlit/secrets.toml", "r", encoding="utf-8") as f:
                config = toml.load(f)
                if "ZHIPU_API_KEY" in config:
                    os.environ["ZHIPU_API_KEY"] = config["ZHIPU_API_KEY"]
                    print("✅ 已加载本地 API Key")
    except Exception:
        pass

# ================= 测试用例 (精心构造) =================
TEST_CASES = [
    {
        "id": "CASE_1_APPLICATION",
        "desc": "❌ [假阳性测试] B细胞线粒体钙研究 (应用文)",
        "paper": {
            "title": "Calcium signaling in B cells using genetically encoded sensors.",
            "abstract": "Chemical dyes as well as a genetically encoded Ca2+ sensor with a mitochondrial targeting sequence were used to study mitochondrial Ca2+ dynamics in response to various stimuli. We show that mitochondrial Ca2+ uptake is dependent on MCU.",
            "doi": "test/001"
        }
    },
    {
        "id": "CASE_2_DEVELOPMENT",
        "desc": "✅ [真阳性测试] 典型的探针开发",
        "paper": {
            "title": "Development of a high-performance genetically encoded sensor for dopamine.",
            "abstract": "Here we report the development of GRAB_DA3.0. We screened a linker library and performed site-directed mutagenesis to improve sensitivity by 5-fold compared to previous versions.",
            "doi": "test/002"
        }
    },
    {
        "id": "CASE_3_NO_ABSTRACT",
        "desc": "⚠️ [极限测试] 只有标题 (Cell Discovery 那篇)",
        "paper": {
            "title": "A genetically encoded ratiometric indicator for tryptophan.",
            "abstract": "[No Abstract]",
            "doi": "test/003"
        }
    }
]

def run_lab_test():
    print("\n🔬 === AI 逻辑实验室启动 ===\n")
    
    if not os.environ.get("ZHIPU_API_KEY"):
        print("❌ 错误: 未找到 API Key，请设置环境变量或 secrets.toml")
        return

    for case in TEST_CASES:
        print(f"🧪 正在测试: {case['desc']}")
        print(f"   📄 标题: {case['paper']['title']}")
        
        # 调用 AI
        result = analyze_one_paper(case['paper'])
        
        if result:
            print(f"   🤖 AI 判定: {'✅ NEW PROBE' if result.get('is_new') else '❌ REJECTED'}")
            print(f"   🧠 推理过程: {result.get('reasoning', 'No reasoning provided')}")
            if result.get('is_new'):
                print(f"   🏷️  提取名称: {result.get('probe_name')}")
        else:
            print("   ⚠️ 分析失败 (网络或API错误)")
        
        print("-" * 60)

if __name__ == "__main__":
    run_lab_test()