import streamlit as st
import pandas as pd
import json
import os
import datetime
from github import Github # 用于连接 GitHub

# ================= 1. 配置与初始化 =================
st.set_page_config(
    page_title="Auto Sensor Base",
    page_icon="🧬",
    layout="wide"
)

# 读取 Secrets
GITHUB_TOKEN = st.secrets.get("GITHUB_TOKEN")
REPO_NAME = st.secrets.get("GITHUB_REPO")
ADMIN_PWD = st.secrets.get("ADMIN_PASSWORD")
DATA_FILE = "processed_probes.json"

# ================= 2. GitHub 同步功能 (核心) =================
def update_github_data(new_data_list):
    """将修改后的数据推送到 GitHub"""
    if not GITHUB_TOKEN or not REPO_NAME:
        st.error("❌ GitHub Token or Repo name not set in Secrets!")
        return False
        
    try:
        g = Github(GITHUB_TOKEN)
        repo = g.get_repo(REPO_NAME)
        contents = repo.get_contents(DATA_FILE)
        
        # 将数据转换为 JSON 字符串
        json_content = json.dumps(new_data_list, indent=4, ensure_ascii=False)
        
        # 提交更新
        repo.update_file(
            path=contents.path,
            message="🤖 Admin: Manual data update via Streamlit",
            content=json_content,
            sha=contents.sha
        )
        return True
    except Exception as e:
        st.error(f"GitHub Sync Error: {e}")
        return False

# ================= 3. 主题配色 =================
def get_theme_config(is_light_mode):
    # (保持你之前的配色代码不变，这里省略以节省篇幅，请直接用你现在的 get_theme_config 函数)
    # ... 请保留你之前那个完美的 get_theme_config 函数 ...
    if is_light_mode:
        return {
            "mode": "light", "bg_color": "#F7F9FB", "sidebar_bg": "#FFFFFF", "text_main": "#2D3748",
            "text_sidebar": "#1A202C", "meta_color": "#718096", "border_color": "#E2E8F0",
            "card_bg": "#FFFFFF", "card_border": "1px solid #E2E8F0", "card_shadow": "0 2px 5px rgba(0,0,0,0.03)",
            "btn_bg": "#FFFFFF", "btn_text": "#4A5568", "btn_border": "#CBD5E0", "btn_hover_bg": "#F7FAFC",
            "badge_target_bg": "#E6FFFA", "badge_target_text": "#2C7A7B", "badge_type_bg": "#EDF2F7", "badge_type_text": "#4A5568",
            "header_visibility": "hidden"
        }
    else:
        return {
            "mode": "dark", "bg_color": "#0E1117", "sidebar_bg": "#262730", "text_main": "#FAFAFA",
            "text_sidebar": "#FAFAFA", "meta_color": "#9CA3AF", "border_color": "#4B5563",
            "card_bg": "#1F2937", "card_border": "1px solid #374151", "card_shadow": "none",
            "btn_bg": "#374151", "btn_text": "#E5E7EB", "btn_border": "#4B5563", "btn_hover_bg": "#4B5563",
            "badge_target_bg": "rgba(16, 185, 129, 0.15)", "badge_target_text": "#6EE7B7", "badge_type_bg": "rgba(255, 255, 255, 0.1)", "badge_type_text": "#D1D5DB",
            "header_visibility": "visible"
        }

def inject_custom_css(t):
    # (保持你之前的 CSS 代码不变)
    st.markdown(f"""
    <style>
        .stApp {{ background-color: {t['bg_color']}; color: {t['text_main']}; }}
        .stDecoration {{ display: none !important; }}
        header[data-testid="stHeader"] {{ background-color: transparent !important; visibility: {t['header_visibility']}; }}
        [data-testid="stSidebar"] {{ background-color: {t['sidebar_bg']}; border-right: 1px solid #E2E8F0; }}
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] span, [data-testid="stSidebar"] div, [data-testid="stSidebar"] label {{ color: {t['text_sidebar']} !important; }}
        .block-container {{ padding-top: 1.5rem; padding-bottom: 3rem; }}
        [data-testid="stVerticalBlockBorderWrapper"] > div {{ background-color: {t['card_bg']}; border: {t['card_border']} !important; box-shadow: {t['card_shadow']}; border-radius: 10px; padding: 1.2rem; }}
        .stButton button, [data-testid="stLinkButton"] a {{ background-color: {t['btn_bg']} !important; color: {t['btn_text']} !important; border: 1px solid {t['btn_border']} !important; border-radius: 6px; }}
        .probe-title {{ font-size: 1.2rem; font-weight: 700; margin-bottom: 6px; color: {t['text_main']}; }}
        .probe-meta {{ font-size: 0.9rem; color: {t['meta_color']}; font-family: 'Source Sans Pro', sans-serif; }}
        .badge-target {{ background-color: {t['badge_target_bg']}; color: {t['badge_target_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; font-weight: 600; display: inline-block; margin-right: 6px; }}
        .badge-type {{ background-color: {t['badge_type_bg']}; color: {t['badge_type_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; display: inline-block; }}
        .badge-new {{ background: linear-gradient(135deg, #FFD700 0%, #F59E0B 100%); color: white; padding: 2px 6px; border-radius: 4px; font-size: 0.7rem; font-weight: 800; margin-left: 8px; vertical-align: middle; }}
    </style>
    """, unsafe_allow_html=True)

# ================= 4. 数据加载 =================
def load_data():
    # 优先读 GitHub 也可以，但为了速度，Streamlit 会读本地副本。
    # 每次更新后，Streamlit Cloud 会自动 reload，所以读本地即可。
    if not os.path.exists(DATA_FILE):
        return []
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data # 返回 List 而不是 DataFrame，方便操作

def get_color_circle_html(color_name):
    # (保持不变)
    c = str(color_name).lower()
    hex_color = "#D1D5DB"
    if "green" in c: hex_color = "#10B981"
    elif "red" in c: hex_color = "#EF4444"
    elif "blue" in c or "cyan" in c: hex_color = "#3B82F6"
    elif "yellow" in c or "gold" in c: hex_color = "#F59E0B"
    elif "orange" in c: hex_color = "#F97316"
    elif "purple" in c: hex_color = "#8B5CF6"
    return f"""<div style="width: 14px; height: 14px; background-color: {hex_color}; border-radius: 50%; display: inline-block; box-shadow: 0 0 0 2px {hex_color}30;"></div>"""

# ================= 5. 管理员功能模块 =================

def render_admin_panel(current_data):
    """侧边栏：添加新探针"""
    with st.sidebar:
        st.markdown("---")
        st.markdown("### 🛠️ Admin Panel")
        
        with st.expander("➕ Manual Add Probe", expanded=False):
            with st.form("add_probe_form"):
                new_name = st.text_input("Probe Name", placeholder="e.g. Mel-G")
                new_target = st.text_input("Target", placeholder="e.g. Melatonin")
                new_color = st.selectbox("Color", ["Green", "Red", "Blue", "Yellow", "Orange", "Other"])
                new_type = st.text_input("Type", placeholder="e.g. GPCR-based")
                new_title = st.text_input("Paper Title")
                new_journal = st.text_input("Journal")
                new_year = st.text_input("Year", value=str(datetime.datetime.now().year))
                new_doi = st.text_input("DOI Link", placeholder="https://doi.org/...")
                new_abstract = st.text_area("Abstract")
                
                submitted = st.form_submit_button("Submit & Save")
                
                if submitted:
                    if not new_name or not new_title:
                        st.error("Name and Title are required!")
                    else:
                        new_entry = {
                            "probe_name": new_name,
                            "target": new_target,
                            "color": new_color,
                            "type": new_type,
                            "title": new_title,
                            "journal": new_journal,
                            "date": new_year,
                            "doi": new_doi,
                            "abstract": new_abstract,
                            "is_new": True # 手动添加的默认是新
                        }
                        
                        # 插入到最前面
                        current_data.insert(0, new_entry)
                        
                        with st.spinner("Syncing to GitHub..."):
                            if update_github_data(current_data):
                                st.success("Added successfully!")
                                st.rerun() # 刷新页面
                            else:
                                st.error("Failed to sync.")

# ================= 6. 渲染主列表 (带删除功能) =================

def render_main_feed(data, theme, is_admin):
    st.header("🚀 Latest Probes")

    if not data:
        st.info("No data available.")
        return

    # 转换为 DataFrame 方便排序，但在操作删除时我们要用原始 list
    # 这里为了简单，我们直接遍历 List
    # 如果需要排序，建议先对 list 排序
    data.sort(key=lambda x: str(x.get('date', '0')), reverse=True)

    for index, row in enumerate(data):
        # 布局
        with st.container(border=True):
            c1, c2, c3 = st.columns([0.2, 5.5, 1.0])
            
            with c1:
                st.markdown(f"<div style='padding-top: 4px; text-align: center;'>{get_color_circle_html(row.get('color', ''))}</div>", unsafe_allow_html=True)
            
            with c2:
                # 标题 + 徽章
                current_year = datetime.datetime.now().year
                pub_year = str(row.get('date', ''))
                is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
                new_badge = '<span class="badge-new">NEW</span>' if is_new else ""
                
                st.markdown(f"""
                <div class="probe-title">{row.get('probe_name', 'Unknown')} {new_badge}</div>
                <div style="margin-top: 8px; line-height: 1.8;">
                    <span class="badge-target">{row.get('target', 'N/A')}</span>
                    <span class="badge-type">{row.get('type', 'N/A')}</span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{row.get('journal', 'Unknown')}</i></span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {row.get('date', 'N/A')}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3:
                # 这里的 UI 逻辑：
                # 如果是管理员，显示 "删除" 按钮
                # 如果是普通用户，显示 "Read" 按钮
                if is_admin:
                    if st.button("🗑️ Delete", key=f"del_{index}", type="primary", use_container_width=True):
                        # 删除逻辑
                        data.pop(index)
                        with st.spinner("Deleting & Syncing..."):
                            if update_github_data(data):
                                st.success("Deleted!")
                                st.rerun()
                else:
                    st.markdown("<div style='height: 6px'></div>", unsafe_allow_html=True) 
                    if row.get('doi') and "http" in row['doi']:
                        st.link_button("Read", row['doi'], use_container_width=True)
                    else:
                        st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # 摘要
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)

# ================= 7. 程序入口 =================
def main():
    # 1. 认证逻辑
    # 在 Sidebar 底部放一个密码框，如果输入正确，开启 Admin 模式
    # 先渲染 Sidebar 上半部分
    with st.sidebar:
        st.title("🧬 Auto Sensor Base")
        st.caption("Automated Tracking System")
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
        
    theme_config = get_theme_config(is_light)
    inject_custom_css(theme_config)
    
    # 加载数据
    data_list = load_data()
    
    # Sidebar 内容
    with st.sidebar:
        st.markdown("<br>", unsafe_allow_html=True)
        if data_list:
            df_temp = pd.DataFrame(data_list)
            csv = df_temp.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Download CSV", csv, "probes.csv", "text/csv", use_container_width=True)
        st.divider()
        
        # 筛选器逻辑 (略微简化，逻辑同前)
        # ... (如果你需要筛选器，把之前的代码贴回这里) ...
        
        # === 🔐 管理员认证入口 ===
        st.markdown("<br><br>", unsafe_allow_html=True)
        with st.expander("🔐 Admin Login"):
            pwd = st.text_input("Password", type="password")
            is_admin = (pwd == ADMIN_PWD)
            if is_admin:
                st.success("Admin Mode Active")
            elif pwd:
                st.error("Wrong Password")
    
    # 如果是管理员，渲染管理面板
    if is_admin:
        render_admin_panel(data_list)
        
    # 渲染主界面
    render_main_feed(data_list, theme_config, is_admin)

if __name__ == "__main__":
    main()