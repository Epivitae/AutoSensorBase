import streamlit as st
import pandas as pd
import json
import os
import datetime
from github import Github

# ================= 1. 配置与初始化 =================
st.set_page_config(
    page_title="Auto Sensor Base",
    page_icon="🧬",
    layout="wide"
)

# 读取 Secrets (云端配置)
GITHUB_TOKEN = st.secrets.get("GITHUB_TOKEN")
REPO_NAME = st.secrets.get("GITHUB_REPO")
ADMIN_PWD = st.secrets.get("ADMIN_PASSWORD")
DATA_FILE = "processed_probes.json"

# ================= 2. GitHub 同步功能 =================
def update_github_data(new_data_list):
    """将修改后的数据推送到 GitHub"""
    if not GITHUB_TOKEN or not REPO_NAME:
        st.error("❌ GitHub Token or Repo name not set in Secrets!")
        return False
        
    try:
        g = Github(GITHUB_TOKEN)
        repo = g.get_repo(REPO_NAME)
        contents = repo.get_contents(DATA_FILE)
        
        json_content = json.dumps(new_data_list, indent=4, ensure_ascii=False)
        
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

# ================= 3. 主题配色定义 =================
def get_theme_config(is_light_mode):
    if is_light_mode:
        # --- 🌞 白天模式 ---
        return {
            "mode": "light", "bg_color": "#F7F9FB", "sidebar_bg": "#FFFFFF", "text_main": "#2D3748",
            "text_sidebar": "#1A202C", "meta_color": "#718096", "border_color": "#E2E8F0",
            "card_bg": "#FFFFFF", "card_border": "1px solid #E2E8F0", "card_shadow": "0 2px 5px rgba(0,0,0,0.03)",
            "btn_bg": "#FFFFFF", "btn_text": "#4A5568", "btn_border": "#CBD5E0", "btn_hover_bg": "#F7FAFC",
            "badge_target_bg": "#E6FFFA", "badge_target_text": "#2C7A7B", "badge_type_bg": "#EDF2F7", "badge_type_text": "#4A5568",
            "header_visibility": "hidden"
        }
    else:
        # --- 🌜 黑夜模式 ---
        return {
            "mode": "dark", "bg_color": "#0E1117", "sidebar_bg": "#262730", "text_main": "#FAFAFA",
            "text_sidebar": "#FAFAFA", "meta_color": "#9CA3AF", "border_color": "#4B5563",
            "card_bg": "#1F2937", "card_border": "1px solid #374151", "card_shadow": "none",
            "btn_bg": "#374151", "btn_text": "#E5E7EB", "btn_border": "#4B5563", "btn_hover_bg": "#4B5563",
            "badge_target_bg": "rgba(16, 185, 129, 0.15)", "badge_target_text": "#6EE7B7", "badge_type_bg": "rgba(255, 255, 255, 0.1)", "badge_type_text": "#D1D5DB",
            "header_visibility": "visible"
        }

# ================= 4. 样式注入 =================
def inject_custom_css(t):
    st.markdown(f"""
    <style>
        .stApp {{ background-color: {t['bg_color']}; color: {t['text_main']}; }}
        .stDecoration {{ display: none !important; }}
        header[data-testid="stHeader"] {{ background-color: transparent !important; visibility: {t['header_visibility']}; }}
        
        [data-testid="stSidebar"] {{ background-color: {t['sidebar_bg']}; border-right: 1px solid #E2E8F0; }}
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] span, 
        [data-testid="stSidebar"] div, [data-testid="stSidebar"] label, [data-testid="stSidebar"] p {{ color: {t['text_sidebar']} !important; }}
        [data-testid="stSidebar"] a {{ color: {t['text_sidebar']} !important; opacity: 0.8; }}
        [data-testid="stSidebar"] .stSelectbox > div > div {{ background-color: {t['bg_color']}; color: {t['text_main']}; border-color: {t['btn_border']}; }}

        .block-container {{ padding-top: 1.5rem; padding-bottom: 3rem; }}

        [data-testid="stVerticalBlockBorderWrapper"] > div {{
            background-color: {t['card_bg']}; border: {t['card_border']} !important;
            box-shadow: {t['card_shadow']}; border-radius: 10px; padding: 1.2rem;
        }}
        
        .stButton button, [data-testid="stLinkButton"] a, [data-testid="stDownloadButton"] button {{
            background-color: {t['btn_bg']} !important; color: {t['btn_text']} !important;
            border: 1px solid {t['btn_border']} !important; border-radius: 6px; font-weight: 500;
            transition: all 0.2s ease; box-shadow: 0 1px 2px rgba(0,0,0,0.05); text-decoration: none !important;
        }}
        .stButton button:hover, [data-testid="stLinkButton"] a:hover, [data-testid="stDownloadButton"] button:hover {{
            background-color: {t['btn_hover_bg']} !important; border-color: #A0AEC0 !important;
            transform: translateY(-1px); color: {t['text_main']} !important;
        }}

        .probe-title {{ font-size: 1.2rem; font-weight: 700; margin-bottom: 6px; color: {t['text_main']}; letter-spacing: -0.01em; }}
        .probe-meta {{ font-size: 0.9rem; color: {t['meta_color']}; font-family: 'Source Sans Pro', sans-serif; }}

        .badge-target {{ background-color: {t['badge_target_bg']}; color: {t['badge_target_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; font-weight: 600; display: inline-block; margin-right: 6px; }}
        .badge-type {{ background-color: {t['badge_type_bg']}; color: {t['badge_type_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; display: inline-block; }}
        .badge-new {{ background: linear-gradient(135deg, #FFD700 0%, #F59E0B 100%); color: white; padding: 2px 6px; border-radius: 4px; font-size: 0.7rem; font-weight: 800; margin-left: 8px; vertical-align: middle; box-shadow: 0 2px 5px rgba(245, 158, 11, 0.4); }}
    </style>
    """, unsafe_allow_html=True)

# ================= 5. 数据处理 =================
def load_data():
    if not os.path.exists(DATA_FILE):
        return []
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data

def get_color_circle_html(color_name):
    c = str(color_name).lower()
    hex_color = "#D1D5DB"
    if "green" in c: hex_color = "#10B981"
    elif "red" in c: hex_color = "#EF4444"
    elif "blue" in c or "cyan" in c: hex_color = "#3B82F6"
    elif "yellow" in c or "gold" in c: hex_color = "#F59E0B"
    elif "orange" in c: hex_color = "#F97316"
    elif "purple" in c: hex_color = "#8B5CF6"
    return f"""<div style="width: 14px; height: 14px; background-color: {hex_color}; border-radius: 50%; display: inline-block; box-shadow: 0 0 0 2px {hex_color}30;"></div>"""

def extract_years(data):
    if not data: return "N/A", "N/A"
    try:
        df = pd.DataFrame(data)
        years = df['date'].astype(str).str.extract(r'(\d{4})').astype(float)
        return int(years.min().iloc[0]), int(years.max().iloc[0])
    except:
        return "2021", datetime.datetime.now().year

# ================= 6. 管理员面板 =================
def render_admin_panel(current_data):
    with st.sidebar:
        st.markdown("---")
        st.markdown("### 🛠️ Admin Panel")
        with st.expander("➕ Manual Add Probe", expanded=True):
            with st.form("add_probe_form"):
                new_name = st.text_input("Name")
                new_target = st.text_input("Target")
                new_color = st.selectbox("Color", ["Green", "Red", "Blue", "Yellow", "Orange"])
                new_doi = st.text_input("DOI/Link")
                new_year = st.text_input("Year", value=str(datetime.datetime.now().year))
                submitted = st.form_submit_button("Submit")
                if submitted:
                    new_entry = {
                        "probe_name": new_name, "target": new_target, "color": new_color,
                        "doi": new_doi, "date": new_year, "is_new": True,
                        "type": "Manual", "journal": "Manual Entry", "abstract": "Manually added."
                    }
                    current_data.insert(0, new_entry)
                    if update_github_data(current_data):
                        st.success("Added!"); st.rerun()

# ================= 7. 渲染组件 =================

def render_sidebar(data, theme):
    with st.sidebar:
        st.title("Auto Sensor Base")
        st.caption("Automated Tracking System")
        
        # 模式切换
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
        st.markdown("<br>", unsafe_allow_html=True)

        if data:
            df = pd.DataFrame(data)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Download CSV Dataset", csv, "probes.csv", "text/csv", use_container_width=True)
        st.divider()

        filtered = list(data)
        if data:
            df = pd.DataFrame(data)
            all_targets = ["All"] + sorted(list(df['target'].astype(str).unique()))
            all_colors = ["All"] + sorted(list(df['color'].astype(str).unique()))
            
            sel_target = st.selectbox("Target Molecule", all_targets)
            sel_color = st.selectbox("Fluorescence Color", all_colors)
            
            if sel_target != "All": filtered = [d for d in filtered if str(d.get('target')) == sel_target]
            if sel_color != "All": filtered = [d for d in filtered if str(d.get('color')) == sel_color]

        # 统计面板
        min_y, max_y = extract_years(filtered)
        st.markdown(f"""
        <div style='margin-top: 20px; padding: 15px; background: rgba(0,0,0,0.03); border-radius: 12px; text-align: center; border: 1px solid {theme['border_color']};'>
            <div style='font-size: 0.8rem; color: {theme['meta_color']}; letter-spacing: 1px;'>TOTAL SENSORS</div>
            <div style='font-size: 2rem; font-weight: 800; color: {theme['text_sidebar']};'>{len(filtered)}</div>
            <div style='font-size: 0.8rem; color: {theme['text_sidebar']}; opacity: 0.8;'>({min_y} - {max_y})</div>
        </div>
        """, unsafe_allow_html=True)
        
        # === 🔐 管理员登录 (最关键的部分) ===
        st.markdown("<br>", unsafe_allow_html=True)
        is_admin = False
        with st.expander("🔐 Admin Login"):
            pwd = st.text_input("Password", type="password")
            if pwd == ADMIN_PWD and ADMIN_PWD:
                is_admin = True
                st.success("Mode: Admin")
        
        # 底部品牌 Footer
        st.markdown(f"""
        <div style='margin-top: 40px; padding-top: 20px; border-top: 1px solid {theme['border_color']}; text-align: center;'>
            <div style='font-weight: 600; font-size: 0.9rem; margin-bottom: 4px; color: {theme['text_sidebar']};'>Chimera Nano Sensor Team</div>
            <div style='font-size: 0.8rem; opacity: 0.8; margin-bottom: 8px; color: {theme['text_sidebar']};'>Institute of Neuroscience, CAS</div>
            <div style='font-size: 0.8rem;'>
                <a href="http://www.cns.ac.cn" target="_blank" style='text-decoration: none; border-bottom: 1px dotted; color: {theme['text_sidebar']};'>www.cns.ac.cn</a>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        return is_light, filtered, is_admin

def render_main_feed(data, theme, is_admin):
    st.header("🚀 Latest Probes")

    if not data:
        st.info("No data available.")
        return

    # 简单排序
    data.sort(key=lambda x: str(x.get('date', '0')), reverse=True)

    for index, row in enumerate(data):
        current_year = datetime.datetime.now().year
        pub_year = str(row.get('date', ''))
        is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
        new_badge = '<span class="badge-new">NEW</span>' if is_new else ""

        with st.container(border=True):
            c1, c2, c3 = st.columns([0.2, 5.5, 1.0])
            with c1: st.markdown(f"<div style='padding-top: 4px; text-align: center;'>{get_color_circle_html(row.get('color',''))}</div>", unsafe_allow_html=True)
            with c2:
                st.markdown(f"""
                <div class="probe-title">{row.get('probe_name','Unknown')} {new_badge}</div>
                <div style="margin-top: 8px; line-height: 1.8;">
                    <span class="badge-target">{row.get('target','N/A')}</span>
                    <span class="badge-type">{row.get('type','N/A')}</span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{row.get('journal','Unknown')}</i></span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {row.get('date','N/A')}</span>
                </div>
                """, unsafe_allow_html=True)
            with c3:
                st.markdown("<div style='height: 6px'></div>", unsafe_allow_html=True)
                
                # === 按钮逻辑：管理员显示删除，普通人显示阅读 ===
                if is_admin:
                    if st.button("🗑️ Delete", key=f"del_{index}", type="primary", use_container_width=True):
                        data.pop(index) # 删除该条
                        if update_github_data(data):
                            st.success("Deleted!"); st.rerun()
                else:
                    if row.get('doi') and "http" in row['doi']:
                        st.link_button("Read", row['doi'], use_container_width=True)
                    else:
                        st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)

# ================= 8. 程序入口 =================
def main():
    # 1. 加载数据
    data_list = load_data()
    
    # 2. 渲染 Sidebar 并获取状态 (包括是否管理员)
    # 我们先临时获取 theme 来渲染 sidebar
    temp_theme = get_theme_config(False) 
    
    # 重新组织：为了让 toggle 决定 theme，我们需要先渲染 toggle
    # 但 sidebar 的其他内容需要 theme。
    # 解决方案：分两步渲染 Sidebar
    
    with st.sidebar:
        # Step 1: 只有 Toggle
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
    
    # 获取真实 Theme
    theme_config = get_theme_config(is_light)
    inject_custom_css(theme_config)
    
    # Step 2: 渲染 Sidebar 剩余部分 (Filters, Stats, Admin, Footer)
    _, filtered_data, is_admin = render_sidebar(data_list, theme_config)
    
    # 3. 如果是管理员，额外显示添加面板
    if is_admin:
        render_admin_panel(data_list)
    
    # 4. 渲染主列表
    render_main_feed(filtered_data, theme_config, is_admin)

if __name__ == "__main__":
    main()