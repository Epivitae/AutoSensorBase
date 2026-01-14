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
# 在侧边栏顶部添加切换开关
with st.sidebar:
    st.title("🧬 Auto-DB")
    # 默认为 False (即黑夜模式，符合科研习惯)
    is_light_mode = st.toggle("🌞 Day Mode / 🌜 Night", value=False)

# 定义两套配色方案
if is_light_mode:
    # --- 白天模式配色 ---
    theme = {
        "bg_color": "#ffffff",
        "sidebar_bg": "#f8f9fa",
        "text_color": "#31333F",
        "card_bg": "#ffffff",
        "card_border": "#e0e0e0",
        "card_shadow": "0 2px 5px rgba(0,0,0,0.05)",
        "badge_target_bg": "#e6fffa",
        "badge_target_text": "#047857",
        "badge_target_border": "#b7ebd6",
        "badge_type_bg": "#f1f3f5",
        "badge_type_text": "#495057",
        "badge_type_border": "#dee2e6",
        "meta_color": "#666666"
    }
else:
    # --- 黑夜模式配色 (默认) ---
    theme = {
        "bg_color": "#0e1117",
        "sidebar_bg": "#262730",
        "text_color": "#fafafa",
        "card_bg": "#1f2937", #稍微比背景亮一点
        "card_border": "#374151",
        "card_shadow": "none",
        "badge_target_bg": "rgba(16, 185, 129, 0.2)", # 半透明绿
        "badge_target_text": "#6ee7b7", # 荧光绿
        "badge_target_border": "transparent",
        "badge_type_bg": "rgba(255, 255, 255, 0.1)",
        "badge_type_text": "#d1d5db",
        "badge_type_border": "transparent",
        "meta_color": "#9ca3af"
    }

# ================= CSS Injection =================
# 将 Python 变量注入到 CSS 中
st.markdown(f"""
<style>
    /* 全局背景和文字 */
    .stApp {{
        background-color: {theme['bg_color']};
        color: {theme['text_color']};
    }}
    
    /* 侧边栏背景 */
    [data-testid="stSidebar"] {{
        background-color: {theme['sidebar_bg']};
    }}

    /* 调整顶部留白 */
    .block-container {{
        padding-top: 2rem;
        padding-bottom: 3rem;
    }}

    /* 卡片容器样式 (Streamlit 的 st.container(border=True) 对应的类) */
    [data-testid="stVerticalBlockBorderWrapper"] > div {{
        background-color: {theme['card_bg']};
        border-color: {theme['card_border']} !important;
        box-shadow: {theme['card_shadow']};
        transition: transform 0.2s ease, box-shadow 0.2s ease;
    }}
    /* 卡片悬停效果（仅白天模式明显） */
    [data-testid="stVerticalBlockBorderWrapper"] > div:hover {{
        border-color: #888 !important;
    }}

    /* 标题样式 */
    .probe-title {{
        font-size: 1.15rem;
        font-weight: 700;
        margin-bottom: 2px;
        color: {theme['text_color']};
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    }}

    /* 元数据文本 */
    .probe-meta {{
        font-size: 0.85rem;
        color: {theme['meta_color']};
        font-weight: 400;
    }}

    /* 徽章：Target */
    .badge-target {{
        background-color: {theme['badge_target_bg']};
        color: {theme['badge_target_text']};
        padding: 2px 8px;
        border-radius: 6px;
        font-size: 0.75rem;
        font-weight: 600;
        border: 1px solid {theme['badge_target_border']};
        display: inline-block;
        margin-right: 6px;
    }}

    /* 徽章：Type */
    .badge-type {{
        background-color: {theme['badge_type_bg']};
        color: {theme['badge_type_text']};
        padding: 2px 8px;
        border-radius: 6px;
        font-size: 0.75rem;
        border: 1px solid {theme['badge_type_border']};
        display: inline-block;
    }}

    /* NEW 星标徽章 (保持金黄色) */
    .badge-new {{
        background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);
        color: #fff;
        padding: 2px 6px;
        border-radius: 4px;
        font-size: 0.65rem;
        font-weight: bold;
        text-shadow: 0 1px 1px rgba(0,0,0,0.2);
        margin-left: 8px;
        vertical-align: middle;
        letter-spacing: 0.5px;
    }}

    /* 按钮紧凑化 */
    .stButton button {{
        height: 2.2rem;
        line-height: 1;
        border-radius: 8px;
        border: 1px solid {theme['card_border']};
    }}
    
    /* 去除链接下划线 */
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
    """Returns a colored SVG circle with a subtle shadow"""
    c = str(color_name).lower()
    hex_color = "#ccc" 
    if "green" in c: hex_color = "#10b981" # Emerald 500
    elif "red" in c: hex_color = "#ef4444" # Red 500
    elif "blue" in c or "cyan" in c: hex_color = "#3b82f6" # Blue 500
    elif "yellow" in c or "gold" in c: hex_color = "#eab308" # Yellow 500
    elif "orange" in c: hex_color = "#f97316" # Orange 500
    elif "purple" in c: hex_color = "#a855f7" # Purple 500
    
    return f"""
    <div style="
        width: 14px; 
        height: 14px; 
        background-color: {hex_color}; 
        border-radius: 50%; 
        display: inline-block;
        box-shadow: 0 0 8px {hex_color}80; /* Glow effect */
        vertical-align: middle;
        margin-top: 2px;
    "></div>
    """

# Load Data
df = load_data()

# ================= Sidebar Content =================
with st.sidebar:
    st.caption("Automated Tracking System")
    
    # Download Button
    if not df.empty:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Dataset (CSV)",
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
        
        # Sort filters
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
        
        st.markdown(f"<br><div style='text-align: center; color: {theme['meta_color']}'>Showing <b>{len(filtered_df)}</b> probes</div>", unsafe_allow_html=True)
    else:
        filtered_df = pd.DataFrame()

# ================= Main Content =================
st.header("🚀 Latest Probes")

if filtered_df.empty:
    st.info("No data available yet. Please run the backend script.")
else:
    # Sort by date
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
            c1, c2, c3 = st.columns([0.25, 6, 1.2])
            
            with c1:
                # Color Indicator
                st.markdown(f"<div style='padding-top: 6px;'>{color_circle(row['color'])}</div>", unsafe_allow_html=True)
            
            with c2:
                # Title
                st.markdown(f"""
                <div class="probe-title">
                    {row['probe_name']} {new_badge}
                </div>
                """, unsafe_allow_html=True)
                
                # Metadata line (Compact)
                target = row['target']
                ptype = row.get('type', 'Unknown')
                journal = row.get('journal', 'Unknown Journal')
                date = row.get('date', 'N/A')
                
                st.markdown(f"""
                <div style="margin-top: 6px; line-height: 1.6;">
                    <span class="badge-target">{target}</span>
                    <span class="badge-type">{ptype}</span>
                    <span style="color: {theme['card_border']}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{journal}</i></span>
                    <span style="color: {theme['card_border']}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {date}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3:
                # Button
                st.markdown("<div style='height: 4px'></div>", unsafe_allow_html=True) # Spacer
                if row.get('doi') and "http" in row['doi']:
                    st.link_button("Read Paper", row['doi'], use_container_width=True)
                else:
                    st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # Abstract (Optional)
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='font-size: 0.9rem; color: {theme['text_color']}; opacity: 0.9;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)