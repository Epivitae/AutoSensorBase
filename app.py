import streamlit as st
import pandas as pd
import json
import os
import datetime

# ================= Page Config =================
st.set_page_config(
    page_title="FP-Sensor Auto-DB",
    page_icon="🧬",
    layout="wide"
)

# ================= Theme Logic (核心配色逻辑) =================
with st.sidebar:
    st.title("🧬 Auto-DB")
    # 开关：False=Dark(默认), True=Light
    is_light_mode = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)

# 定义配色方案
if is_light_mode:
    # --- 🌞 白天模式 (Light Mode) ---
    theme = {
        "bg_color": "#F7F9FB",         # 极浅的灰蓝色背景，护眼且高级
        "sidebar_bg": "#FFFFFF",       # 纯白侧边栏
        "sidebar_text": "#1A202C",     # 侧边栏文字深黑
        "main_text": "#2D3748",        # 主内容文字深灰
        "card_bg": "#FFFFFF",          # 卡片纯白
        "card_border": "1px solid #E2E8F0", # 极浅的边框
        "card_shadow": "0 2px 4px rgba(0,0,0,0.02), 0 1px 2px rgba(0,0,0,0.03)", # 几乎不可见的微阴影
        "header_visibility": "hidden", # 隐藏原来的深色 Header
        "badge_target_bg": "#E6FFFA",  # 清新薄荷绿
        "badge_target_text": "#2C7A7B",
        "badge_type_bg": "#EDF2F7",    # 浅灰胶囊
        "badge_type_text": "#4A5568",
        "meta_color": "#718096",       # 元数据灰色
        
        # 按钮样式 (轻量化)
        "btn_bg": "#FFFFFF",
        "btn_text": "#4A5568",
        "btn_border": "#CBD5E0",
        "btn_hover_bg": "#F7FAFC",
        "btn_hover_border": "#A0AEC0"
    }
else:
    # --- 🌜 黑夜模式 (Dark Mode) ---
    theme = {
        "bg_color": "#0E1117",
        "sidebar_bg": "#262730",
        "sidebar_text": "#FAFAFA",
        "main_text": "#FAFAFA",
        "card_bg": "#1F2937",
        "card_border": "1px solid #374151",
        "card_shadow": "none",
        "header_visibility": "visible",
        "badge_target_bg": "rgba(16, 185, 129, 0.15)",
        "badge_target_text": "#6EE7B7",
        "badge_type_bg": "rgba(255, 255, 255, 0.1)",
        "badge_type_text": "#D1D5DB",
        "meta_color": "#9CA3AF",
        
        # 按钮样式 (深色)
        "btn_bg": "#374151",
        "btn_text": "#E5E7EB",
        "btn_border": "#4B5563",
        "btn_hover_bg": "#4B5563",
        "btn_hover_border": "#9CA3AF"
    }

# ================= CSS Injection (黑魔法区域) =================
st.markdown(f"""
<style>
    /* 1. 全局背景与文字 */
    .stApp {{
        background-color: {theme['bg_color']};
        color: {theme['main_text']};
    }}

    /* 2. 隐藏 Streamlit 顶部的彩虹条和深色 Header (针对白天模式优化) */
    header[data-testid="stHeader"] {{
        background-color: transparent !important;
        visibility: {theme.get('header_visibility', 'visible')};
    }}
    /* 彻底隐藏顶部的彩虹装饰条 */
    .stDecoration {{
        display: none !important;
    }}

    /* 3. 侧边栏样式 (强制覆盖) */
    [data-testid="stSidebar"] {{
        background-color: {theme['sidebar_bg']};
        border-right: 1px solid #E2E8F0;
    }}
    /* 强制侧边栏内所有元素（标题、文本、Label）的颜色 */
    [data-testid="stSidebar"] h1, 
    [data-testid="stSidebar"] h2, 
    [data-testid="stSidebar"] h3, 
    [data-testid="stSidebar"] span, 
    [data-testid="stSidebar"] label,
    [data-testid="stSidebar"] div {{
        color: {theme['sidebar_text']} !important;
    }}
    /* 侧边栏下拉框的优化 */
    [data-testid="stSidebar"] [data-baseweb="select"] div {{
        background-color: {theme['bg_color']};
        color: {theme['main_text']};
    }}

    /* 4. 顶部留白调整 (去除 Header 后的补位) */
    .block-container {{
        padding-top: 1rem;
        padding-bottom: 3rem;
    }}

    /* 5. 卡片样式 */
    [data-testid="stVerticalBlockBorderWrapper"] > div {{
        background-color: {theme['card_bg']};
        border: {theme['card_border']} !important;
        box-shadow: {theme['card_shadow']};
        border-radius: 10px;
        padding: 1.2rem;
    }}

    /* 6. 按钮样式 (Read 按钮) */
    .stButton button, [data-testid="stLinkButton"] {{
        background-color: {theme['btn_bg']} !important;
        color: {theme['btn_text']} !important;
        border: 1px solid {theme['btn_border']} !important;
        border-radius: 6px;
        font-weight: 500;
        transition: all 0.2s ease;
        box-shadow: 0 1px 2px rgba(0,0,0,0.05);
    }}
    
    /* 按钮悬停态 */
    .stButton button:hover, [data-testid="stLinkButton"]:hover {{
        background-color: {theme['btn_hover_bg']} !important;
        border-color: {theme['btn_hover_border']} !important;
        transform: translateY(-1px);
        box-shadow: 0 4px 6px rgba(0,0,0,0.08);
    }}
    
    /* 链接文字颜色修正 */
    [data-testid="stLinkButton"] a {{
        color: {theme['btn_text']} !important;
    }}

    /* 7. 字体与排版优化 */
    .probe-title {{
        font-size: 1.25rem;
        font-weight: 700;
        margin-bottom: 6px;
        color: {theme['main_text']};
        letter-spacing: -0.01em;
    }}
    
    .probe-meta {{
        font-size: 0.9rem;
        color: {theme['meta_color']};
        font-family: 'Source Sans Pro', sans-serif;
    }}

    /* 8. 徽章样式 */
    .badge-target {{
        background-color: {theme['badge_target_bg']};
        color: {theme['badge_target_text']};
        padding: 3px 10px;
        border-radius: 100px;
        font-size: 0.8rem;
        font-weight: 600;
        display: inline-block;
        margin-right: 8px;
    }}
    .badge-type {{
        background-color: {theme['badge_type_bg']};
        color: {theme['badge_type_text']};
        padding: 3px 10px;
        border-radius: 100px;
        font-size: 0.8rem;
        display: inline-block;
    }}
    
    /* NEW 星标 */
    .badge-new {{
        background: linear-gradient(135deg, #FFD700 0%, #F59E0B 100%);
        color: white;
        padding: 2px 6px;
        border-radius: 4px;
        font-size: 0.7rem;
        font-weight: 800;
        margin-left: 8px;
        vertical-align: middle;
        box-shadow: 0 2px 5px rgba(245, 158, 11, 0.4);
    }}
    
    /* 隐藏链接下划线 */
    a {{ text-decoration: none !important; }}
</style>
""", unsafe_allow_html=True)

# File Path
DATA_FILE = "processed_probes.json"

def load_data():
    if not os.path.exists(DATA_FILE):
        return pd.DataFrame()
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    return pd.DataFrame(data)

def color_circle(color_name):
    c = str(color_name).lower()
    hex_color = "#D1D5DB" # Default gray
    
    # 调整过的颜色，让白天模式下也好看
    if "green" in c: hex_color = "#10B981"
    elif "red" in c: hex_color = "#EF4444"
    elif "blue" in c or "cyan" in c: hex_color = "#3B82F6"
    elif "yellow" in c or "gold" in c: hex_color = "#F59E0B"
    elif "orange" in c: hex_color = "#F97316"
    elif "purple" in c: hex_color = "#8B5CF6"
    
    return f"""
    <div style="
        width: 16px; 
        height: 16px; 
        background-color: {hex_color}; 
        border-radius: 50%; 
        display: inline-block;
        box-shadow: 0 0 0 2px {hex_color}40; /* 增加一圈淡色光晕 */
        vertical-align: middle;
        margin-top: 4px;
    "></div>
    """

# Load Data
df = load_data()

# ================= Sidebar =================
with st.sidebar:
    st.caption("Automated Tracking System")
    
    if not df.empty:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download CSV Dataset",
            data=csv,
            file_name='probes_database.csv',
            mime='text/csv',
            use_container_width=True
        )
    
    st.divider()

    # Filters
    if not df.empty:
        df['target'] = df['target'].astype(str)
        df['color'] = df['color'].astype(str)
        df['date'] = df['date'].astype(str)
        
        all_targets = ["All"] + sorted(list(df['target'].unique()))
        selected_target = st.selectbox("Target Molecule", all_targets)
        
        all_colors = ["All"] + sorted(list(df['color'].unique()))
        selected_color = st.selectbox("Fluorescence Color", all_colors)
        
        # Apply Filters
        filtered_df = df.copy()
        if selected_target != "All":
            filtered_df = filtered_df[filtered_df['target'] == selected_target]
        if selected_color != "All":
            filtered_df = filtered_df[filtered_df['color'] == selected_color]
        
        # 侧边栏底部统计
        st.markdown(f"""
        <div style='margin-top: 20px; padding: 10px; background: rgba(0,0,0,0.05); border-radius: 8px; text-align: center;'>
            <div style='font-size: 0.8rem; color: {theme['meta_color']}'>Total Probes</div>
            <div style='font-size: 1.5rem; font-weight: bold; color: {theme['sidebar_text']}'>{len(filtered_df)}</div>
        </div>
        """, unsafe_allow_html=True)
    else:
        filtered_df = pd.DataFrame()

# ================= Main Content =================
st.header("🚀 Latest Probes")

if filtered_df.empty:
    st.info("No data available yet.")
else:
    try:
        filtered_df = filtered_df.sort_values(by='date', ascending=False)
    except:
        pass

    for index, row in filtered_df.iterrows():
        # New Badge Logic
        current_year = datetime.datetime.now().year
        pub_year = str(row.get('date', ''))
        is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
        new_badge = '<span class="badge-new">NEW</span>' if is_new else ""

        # Layout Container
        with st.container(border=True):
            # 列比例优化
            c1, c2, c3 = st.columns([0.2, 6, 1.2])
            
            with c1:
                st.markdown(f"<div style='text-align: center;'>{color_circle(row['color'])}</div>", unsafe_allow_html=True)
            
            with c2:
                # Title
                st.markdown(f"""
                <div class="probe-title">
                    {row['probe_name']} {new_badge}
                </div>
                """, unsafe_allow_html=True)
                
                # Meta Data Row
                target = row['target']
                ptype = row.get('type', 'Unknown')
                journal = row.get('journal', 'Unknown Journal')
                date = row.get('date', 'N/A')
                
                st.markdown(f"""
                <div style="margin-top: 8px; line-height: 1.8;">
                    <span class="badge-target">{target}</span>
                    <span class="badge-type">{ptype}</span>
                    <span style="color: {theme.get('meta_color')}; margin: 0 8px; opacity: 0.5;">|</span>
                    <span class="probe-meta"><i>{journal}</i></span>
                    <span style="color: {theme.get('meta_color')}; margin: 0 8px; opacity: 0.5;">•</span>
                    <span class="probe-meta">📅 {date}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3:
                # Button: 垂直居中
                st.markdown("<div style='height: 8px'></div>", unsafe_allow_html=True) 
                if row.get('doi') and "http" in row['doi']:
                    st.link_button("Read Paper", row['doi'], use_container_width=True)
                else:
                    st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # Abstract Expander
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='font-size: 0.95rem; color: {theme['main_text']}; opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)