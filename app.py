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

# ================= Theme Logic =================
with st.sidebar:
    st.title("🧬 Auto-DB")
    # 开关：False=Dark(默认), True=Light
    is_light_mode = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)

# 定义配色方案
if is_light_mode:
    # --- 🌞 现代清爽白天模式 ---
    theme = {
        "bg_color": "#f3f4f6",         # 浅灰背景，让卡片浮起来
        "sidebar_bg": "#ffffff",       # 纯白侧边栏
        "text_color": "#111827",       # 深灰黑，不刺眼
        "card_bg": "#ffffff",          # 纯白卡片
        "card_border": "transparent",  # 白天模式不需要边框，靠阴影
        "card_shadow": "0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)", # 柔和投影
        "badge_target_bg": "#ecfdf5",  # 清新薄荷绿
        "badge_target_text": "#059669",
        "badge_type_bg": "#f3f4f6",    # 浅灰胶囊
        "badge_type_text": "#4b5563",
        "meta_color": "#6b7280",
        "btn_bg": "#ffffff",           # 按钮白底
        "btn_text": "#374151",         # 按钮灰字
        "btn_border": "#d1d5db",       # 按钮灰边
        "btn_hover": "#f9fafb"
    }
else:
    # --- 🌜 极客深色模式 (保持你喜欢的样子) ---
    theme = {
        "bg_color": "#0e1117",
        "sidebar_bg": "#262730",
        "text_color": "#fafafa",
        "card_bg": "#1f2937",
        "card_border": "#374151",
        "card_shadow": "none",
        "badge_target_bg": "rgba(16, 185, 129, 0.15)",
        "badge_target_text": "#6ee7b7",
        "badge_type_bg": "rgba(255, 255, 255, 0.1)",
        "badge_type_text": "#d1d5db",
        "meta_color": "#9ca3af",
        "btn_bg": "#1f2937",
        "btn_text": "#e5e7eb",
        "btn_border": "#4b5563",
        "btn_hover": "#374151"
    }

# ================= CSS Injection =================
st.markdown(f"""
<style>
    /* 全局背景 */
    .stApp {{
        background-color: {theme['bg_color']};
        color: {theme['text_color']};
    }}
    
    /* 侧边栏 */
    [data-testid="stSidebar"] {{
        background-color: {theme['sidebar_bg']};
        border-right: 1px solid {theme.get('btn_border', 'transparent')};
    }}

    /* 顶部留白 */
    .block-container {{
        padding-top: 2rem;
        padding-bottom: 3rem;
    }}

    /* === 卡片样式优化 === */
    [data-testid="stVerticalBlockBorderWrapper"] > div {{
        background-color: {theme['card_bg']};
        border: 1px solid {theme['card_border']} !important;
        box-shadow: {theme['card_shadow']};
        border-radius: 12px; /* 更圆润的角 */
        padding: 1rem;
        transition: all 0.2s ease;
    }}
    
    /* === 按钮样式重构 (去掉原来的黑砖头) === */
    /* 普通按钮 & 链接按钮 */
    .stButton button, [data-testid="stLinkButton"] {{
        background-color: {theme['btn_bg']} !important;
        color: {theme['btn_text']} !important;
        border: 1px solid {theme['btn_border']} !important;
        border-radius: 8px;
        height: 2.2rem;
        line-height: 2.2rem;
        padding: 0 1rem;
        font-weight: 500;
        transition: all 0.2s;
        box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
    }}
    /* 按钮悬停 */
    .stButton button:hover, [data-testid="stLinkButton"]:hover {{
        background-color: {theme['btn_hover']} !important;
        border-color: #9ca3af !important;
        color: {theme['text_color']} !important;
        transform: translateY(-1px);
    }}
    /* 去除链接按钮的下划线和默认色 */
    [data-testid="stLinkButton"] a {{
        text-decoration: none !important;
        color: {theme['btn_text']} !important;
    }}

    /* === 文本排版 === */
    .probe-title {{
        font-size: 1.2rem;
        font-weight: 700;
        margin-bottom: 4px;
        color: {theme['text_color']};
        letter-spacing: -0.025em;
    }}

    .probe-meta {{
        font-size: 0.85rem;
        color: {theme['meta_color']};
    }}

    /* === 徽章系统 === */
    .badge-target {{
        background-color: {theme['badge_target_bg']};
        color: {theme['badge_target_text']};
        padding: 2px 8px;
        border-radius: 9999px; /* 全圆角胶囊 */
        font-size: 0.75rem;
        font-weight: 600;
        display: inline-block;
        margin-right: 6px;
    }}

    .badge-type {{
        background-color: {theme['badge_type_bg']};
        color: {theme['badge_type_text']};
        padding: 2px 8px;
        border-radius: 9999px;
        font-size: 0.75rem;
        font-weight: 500;
        display: inline-block;
    }}

    /* NEW 星标 (保持金色) */
    .badge-new {{
        background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%);
        color: white;
        padding: 2px 6px;
        border-radius: 4px;
        font-size: 0.65rem;
        font-weight: 800;
        margin-left: 8px;
        vertical-align: middle;
        box-shadow: 0 2px 4px rgba(245, 158, 11, 0.3);
    }}

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
    hex_color = "#d1d5db" 
    if "green" in c: hex_color = "#10b981"
    elif "red" in c: hex_color = "#ef4444"
    elif "blue" in c or "cyan" in c: hex_color = "#3b82f6"
    elif "yellow" in c or "gold" in c: hex_color = "#eab308"
    elif "orange" in c: hex_color = "#f97316"
    elif "purple" in c: hex_color = "#a855f7"
    
    # 给圆点加一点光泽
    return f"""
    <div style="
        width: 16px; 
        height: 16px; 
        background-color: {hex_color}; 
        border-radius: 50%; 
        display: inline-block;
        box-shadow: 0 0 0 2px {hex_color}30; 
        vertical-align: middle;
        margin-top: 2px;
    "></div>
    """

# Load Data
df = load_data()

# ================= Sidebar Content =================
with st.sidebar:
    st.caption("Automated Tracking System")
    
    if not df.empty:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download CSV",
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
        
        st.markdown(f"<br><div style='text-align: center; color: {theme['meta_color']}'>Found <b>{len(filtered_df)}</b> probes</div>", unsafe_allow_html=True)
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

        # Layout
        with st.container(border=True):
            # 调整了列宽比例，让按钮不那么挤
            c1, c2, c3 = st.columns([0.2, 5.5, 1.0])
            
            with c1:
                st.markdown(f"<div style='padding-top: 6px;'>{color_circle(row['color'])}</div>", unsafe_allow_html=True)
            
            with c2:
                # Title
                st.markdown(f"""
                <div class="probe-title">
                    {row['probe_name']} {new_badge}
                </div>
                """, unsafe_allow_html=True)
                
                # Metadata line
                target = row['target']
                ptype = row.get('type', 'Unknown')
                journal = row.get('journal', 'Unknown Journal')
                date = row.get('date', 'N/A')
                
                st.markdown(f"""
                <div style="margin-top: 8px; line-height: 1.6;">
                    <span class="badge-target">{target}</span>
                    <span class="badge-type">{ptype}</span>
                    <span style="color: {theme.get('btn_border')}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{journal}</i></span>
                    <span style="color: {theme.get('btn_border')}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {date}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3:
                # Button
                st.markdown("<div style='height: 4px'></div>", unsafe_allow_html=True) 
                if row.get('doi') and "http" in row['doi']:
                    # 使用 link_button
                    st.link_button("Read", row['doi'], use_container_width=True)
                else:
                    st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # Abstract
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='font-size: 0.9rem; color: {theme['text_color']}; opacity: 0.9; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)