from data_fetcher import fetch_broad_probe_papers
from data_analyzer import analyze_one_paper

def test_fetcher():
    print("🔬 Testing Fetcher Logic...")
    # 只抓 1 天，快速验证
    results = fetch_broad_probe_papers(days_back=1)
    print(f"Found {len(results)} papers.")
    if results:
        print(f"Sample: {results[0]['title']}")
    return results

def test_analyzer(paper_data):
    print("\n🔬 Testing Analyzer Logic...")
    if not paper_data:
        # 造一个假数据测试
        paper_data = {
            "title": "Development of a new red GECI",
            "abstract": "We engineered a novel red fluorescent calcium indicator named R-CaMP9.",
            "doi": "test/123"
        }
    
    result = analyze_one_paper(paper_data)
    print(f"AI Result: {result}")

if __name__ == "__main__":
    # 1. 测试爬虫
    papers = test_fetcher()
    
    # 2. 如果爬到了，就拿第一篇测一下 AI；没爬到就用假数据测 AI
    sample = papers[0] if papers else None
    test_analyzer(sample)