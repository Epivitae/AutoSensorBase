import json
import os
import time
import toml
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TextColumn,
    TimeRemainingColumn,
    MofNCompleteColumn
)
from rich.theme import Theme

# 引入你的业务逻辑模块
from data_fetcher import fetch_broad_probe_papers
from data_analyzer import analyze_one_paper

RAW_FILE = "raw_papers.json"
PROCESSED_FILE = "processed_probes.json"

# === 🟢 初始化 Rich Console ===
custom_theme = Theme({
    "info": "dim cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "probe": "magenta"
})
console = Console(theme=custom_theme, force_terminal=True)

# ================= 🟢 环境变量加载 =================
secrets_path = ".streamlit/secrets.toml"
if os.path.exists(secrets_path):
    try:
        secrets = toml.load(secrets_path)
        if "ZHIPU_API_KEY" in secrets:
            os.environ["ZHIPU_API_KEY"] = secrets["ZHIPU_API_KEY"]
            console.print("[success]✅ Successfully loaded API Key from secrets.toml[/]")
    except Exception as e:
        console.print(f"[error]⚠️ Failed to load secrets: {e}[/]")
# ===================================================

def load_json(filename):
    if os.path.exists(filename):
        with open(filename, "r", encoding="utf-8") as f:
            return json.load(f)
    return []

def save_json(filename, data):
    """保存数据到 JSON 文件"""
    # 为了安全，可以先写到临时文件再重命名（可选），这里保持简单直接写
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)

def normalize_doi(doi_string):
    if not doi_string: return ""
    return doi_string.lower().strip().replace("https://doi.org/", "").replace("http://doi.org/", "")

def print_stats(total, analyzed, pending, batch_size, new_fetched=0):
    """打印漂亮的统计表格"""
    table = Table(title="📊 Analysis Pipeline Status", show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", justify="right", style="green")
    table.add_column("Description", style="dim")

    table.add_row("📚 Total Papers", str(total), "Total documents in raw database")
    table.add_row("🌍 Newly Fetched", str(new_fetched), "Added in this run (PubMed)")
    table.add_row("✅ Already Analyzed", str(analyzed), "Processed in previous runs")
    table.add_row("⏳ Pending Queue", str(pending), "Waiting for AI analysis")
    table.add_row("🚀 Current Batch", str(batch_size), "Will be processed now")

    console.print(Panel(table, expand=False))

def main():
    console.print("[bold blue]🚀 [ETL] Starting Daily Pipeline...[/]")

    # ================= 1. 抓取阶段 =================
    raw_data = load_json(RAW_FILE)
    
    existing_raw_dois = set()
    for item in raw_data:
        existing_raw_dois.add(normalize_doi(item.get('doi')))
    
    console.print("[info]🌍 Fetching from PubMed (5 days)...[/]")
    candidates = fetch_broad_probe_papers(days_back=5)
    
    new_raw_count = 0
    for p in candidates:
        clean_doi = normalize_doi(p['doi'])
        if clean_doi not in existing_raw_dois:
            p['ai_analyzed'] = False 
            p['is_probe'] = False
            raw_data.append(p)
            existing_raw_dois.add(clean_doi)
            new_raw_count += 1
            
    if new_raw_count > 0:
        console.print(f"[success]📥 Staged {new_raw_count} new papers to {RAW_FILE}[/]")
        save_json(RAW_FILE, raw_data)
    else:
        console.print("[dim]💤 No new raw papers found from fetcher.[/]")

    # ================= 2. 统计与准备 =================
    total_docs = len(raw_data)
    analyzed_docs = sum(1 for p in raw_data if p.get('ai_analyzed'))
    pending_list = [p for p in raw_data if not p.get('ai_analyzed')]
    pending_count = len(pending_list)
    
    if not pending_list:
        print_stats(total_docs, analyzed_docs, pending_count, 0, new_raw_count)
        console.print("[bold green]✅ All caught up. Workflow finished.[/]")
        return

    processed_data = load_json(PROCESSED_FILE)
    
    processed_dois_set = set()
    for item in processed_data:
        processed_dois_set.add(normalize_doi(item.get('doi')))
    
    # 设定批处理大小
    BATCH_SIZE = 800 
    batch = pending_list[:BATCH_SIZE]
    
    # === 🔥 展示统计面板 ===
    print_stats(total_docs, analyzed_docs, pending_count, len(batch), new_raw_count)

    analyzed_count = 0
    new_probe_count = 0
    updated_probe_count = 0

    # ================= 3. AI 分析阶段 (Rich Progress) =================
    
    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(), 
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeRemainingColumn(),
        console=console,
        expand=True
    )

    with progress:
        task_id = progress.add_task("[cyan]🤖 AI Analyzing...", total=len(batch))
        
        for paper in batch:
            title_short = paper['title'][:50] + "..." if len(paper['title']) > 50 else paper['title']
            
            try:
                result = analyze_one_paper(paper)
                
                # 更新内存中的状态
                paper['ai_analyzed'] = True
                
                if result and result.get('is_new'):
                    final_entry = {**paper, **result}
                    final_entry.pop('ai_analyzed', None)
                    final_entry.pop('is_probe', None)
                    
                    current_doi = normalize_doi(paper.get('doi'))
                    probe_name = result.get('probe_name', 'Unknown')
                    
                    if current_doi in processed_dois_set:
                        # 🔄 更新
                        console.print(f"  [warning]🔄 UPDATE:[/warning] [bold]{probe_name}[/] (Overwrite)")
                        for idx, item in enumerate(processed_data):
                            if normalize_doi(item.get('doi')) == current_doi:
                                processed_data[idx] = final_entry
                                break
                        updated_probe_count += 1
                    else:
                        # 🎉 新增
                        console.print(f"  [success]🎉 NEW PROBE:[/success] [bold magenta]{probe_name}[/]")
                        console.print(f"     [dim]Title: {title_short}[/]")
                        processed_data.append(final_entry)
                        processed_dois_set.add(current_doi)
                        new_probe_count += 1
                    
                    paper['is_probe'] = True
                else:
                    # console.print(f"  [dim]❌ Rejected: {title_short}[/]") 
                    paper['is_probe'] = False
                
                analyzed_count += 1
                
                # === 💾 关键修改：每 10 个保存一次 Checkpoint ===
                if analyzed_count % 10 == 0:
                    # 更新进度条描述让用户知道正在保存
                    progress.update(task_id, description="[yellow]💾 Saving Checkpoint...[/]")
                    
                    save_json(RAW_FILE, raw_data)
                    save_json(PROCESSED_FILE, processed_data)
                    
                    # 保存完改回原来的描述
                    progress.update(task_id, description="[cyan]🤖 AI Analyzing...[/]")
                
                # 推进进度条
                progress.advance(task_id)
                time.sleep(1) 
                
            except Exception as e:
                console.print(f"  [error]⚠️ Error analyzing {title_short}: {e}[/]")
                progress.advance(task_id)
                continue

    # ================= 4. 最终保存 =================
    # 循环结束后的最后一次保存（防止不是 10 的倍数时丢失最后几条）
    if analyzed_count > 0:
        save_json(RAW_FILE, raw_data)
        save_json(PROCESSED_FILE, processed_data)
        
        summary_table = Table(title="💾 Run Summary", show_header=False, box=None)
        summary_table.add_row("Analyzed", str(analyzed_count), style="blue")
        summary_table.add_row("New Probes Found", str(new_probe_count), style="green bold")
        summary_table.add_row("Existing Updated", str(updated_probe_count), style="yellow")
        
        console.print(Panel(summary_table, title="Pipeline Completed", border_style="green"))

if __name__ == "__main__":
    main()