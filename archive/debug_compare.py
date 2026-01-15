import pandas as pd
from data_fetcher import fetch_recent_papers, fetch_broad_probe_papers

def normalize_doi(doi_url):
    """统一 DOI 格式以便对比 (去掉 http 前缀，转小写)"""
    if not doi_url: return "none"
    # 兼容有些数据可能是 float 类型 (NaN)
    return str(doi_url).lower().replace("https://doi.org/", "").replace("http://doi.org/", "").strip()

def main():
    print("⚔️  开始 A/B 测试：旧爬虫 vs 新爬虫")
    days = 30 

    # 1. 运行旧爬虫 (增加防御性代码：or [])
    print(f"\n👴 正在运行旧版 fetch_recent_papers ({days}天)...")
    old_data = fetch_recent_papers(days_back=days) or [] 
    print(f"   -> 旧版捕获: {len(old_data)} 篇")

    # 2. 运行新爬虫 (增加防御性代码：or [])
    print(f"\n👶 正在运行新版 fetch_broad_probe_papers ({days}天)...")
    new_data = fetch_broad_probe_papers(days_back=days) or []
    print(f"   -> 新版捕获: {len(new_data)} 篇")

    # 3. 数据处理与对比
    old_dois = {normalize_doi(p.get('doi')) for p in old_data}
    
    # 新版可能用大写 DOI 或小写 doi，做个兼容
    new_dois = set()
    for p in new_data:
        d = p.get('DOI') or p.get('doi')
        new_dois.add(normalize_doi(d))

    # 集合运算
    common = old_dois & new_dois         
    only_in_new = new_dois - old_dois    
    only_in_old = old_dois - new_dois    

    # 4. 输出分析报告
    print("\n" + "="*40)
    print("📊  对比分析报告")
    print("="*40)
    print(f"✅ 共同捕获: {len(common)} 篇")
    print(f"🚀 新版新增 (独有): {len(only_in_new)} 篇")
    print(f"⚠️ 旧版独有 (新版漏抓): {len(only_in_old)} 篇")

    # 5. 展示新版多抓到了什么
    if only_in_new:
        print("\n🔎 新版额外发现的高价值论文 (示例前 5 篇):")
        count = 0
        for p in new_data:
            d_norm = normalize_doi(p.get('DOI') or p.get('doi'))
            if d_norm in only_in_new:
                # 打印标题和是否像开发类
                is_dev = p.get('Is_Development_Likely', '')
                # 兼容旧版 title 小写
                title = p.get('Title') or p.get('title')
                print(f"   [{is_dev}] {title}")
                count += 1
                if count >= 5: break
    
    # 6. 检查是否有“倒退”
    if only_in_old:
        print("\n🚨 警告：新版漏掉了旧版能抓到的论文 (需检查原因):")
        for p in old_data:
            d_norm = normalize_doi(p.get('doi'))
            if d_norm in only_in_old:
                title = p.get('Title') or p.get('title')
                print(f"   - {title}")

if __name__ == "__main__":
    main()