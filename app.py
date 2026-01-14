import streamlit as st
import pandas as pd
import json
import os
import datetime
import re

# ================= 1. 配置与初始化 =================
st.set_page_config(
    page_title="Auto Sensor Base",
    page_icon="",
    layout="wide"
)

DATA_FILE = "processed_probes.json"

# ================= 2. 主题配色定义 =================
def get_theme_config(is_light_mode):
    if is_light_mode:
        # --- 🌞 白天模式 ---
        return {
            "mode": "light",
            "bg_color": "#F7F9FB",
            "sidebar_bg": "#FFFFFF",
            "text_main": "#2D3748",
            "text_sidebar": "#1A202C",
            "meta_color": "#718096",
            "border_color": "#E2E8F0",
            
            # 卡片
            "card_bg": "#FFFFFF",
            "card_border": "1px solid #E2E8F0",
            "card_shadow": "0 2px 5px rgba(0,0,0,0.03)",
            
            # 按钮
            "btn_bg": "#FFFFFF",
            "btn_text": "#4A5568",
            "btn_border": "#CBD5E0",
            "btn_hover_bg": "#F7FAFC",
            
            # 徽章
            "badge_target_bg": "#E6FFFA",
            "badge_target_text": "#2C7A7B",
            "badge_type_bg": "#EDF2F7",
            "badge_type_text": "#4A5568",
            
            "header_visibility": "hidden"
        }
    else:
        # --- 🌜 黑夜模式 ---
        return {
            "mode": "dark",
            "bg_color": "#0E1117",
            "sidebar_bg": "#262730",
            "text_main": "#FAFAFA",
            "text_sidebar": "#FAFAFA",
            "meta_color": "#9CA3AF",
            "border_color": "#4B5563",
            
            # 卡片
            "card_bg": "#1F2937",
            "card_border": "1px solid #374151",
            "card_shadow": "none",
            
            # 按钮
            "btn_bg": "#374151",
            "btn_text": "#E5E7EB",
            "btn_border": "#4B5563",
            "btn_hover_bg": "#4B5563",
            
            # 徽章
            "badge_target_bg": "rgba(16, 185, 129, 0.15)",
            "badge_target_text": "#6EE7B7",
            "badge_type_bg": "rgba(255, 255, 255, 0.1)",
            "badge_type_text": "#D1D5DB",
            
            "header_visibility": "visible"
        }

# ================= 3. 样式注入引擎 =================
def inject_custom_css(t):
    st.markdown(f"""
    <style>
        /* 全局背景与文字 */
        .stApp {{ background-color: {t['bg_color']}; color: {t['text_main']}; }}
        
        /* 隐藏装饰条与Header */
        .stDecoration {{ display: none !important; }}
        header[data-testid="stHeader"] {{ background-color: transparent !important; visibility: {t['header_visibility']}; }}
        
        /* 侧边栏样式 */
        [data-testid="stSidebar"] {{ background-color: {t['sidebar_bg']}; border-right: 1px solid #E2E8F0; }}
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3,
        [data-testid="stSidebar"] span, [data-testid="stSidebar"] label, [data-testid="stSidebar"] div, [data-testid="stSidebar"] p {{
            color: {t['text_sidebar']} !important;
        }}
        [data-testid="stSidebar"] a {{ color: {t['text_sidebar']} !important; opacity: 0.8; }}
        
        /* 侧边栏输入框 */
        [data-testid="stSidebar"] .stSelectbox > div > div {{
            background-color: {t['bg_color']}; color: {t['text_main']}; border-color: {t['btn_border']};
        }}

        /* 布局调整 */
        .block-container {{ padding-top: 1.5rem; padding-bottom: 3rem; }}

        /* 卡片样式 */
        [data-testid="stVerticalBlockBorderWrapper"] > div {{
            background-color: {t['card_bg']}; border: {t['card_border']} !important;
            box-shadow: {t['card_shadow']}; border-radius: 10px; padding: 1.2rem;
        }}
        
        /* 按钮通用样式 */
        .stButton button, [data-testid="stLinkButton"] a, [data-testid="stDownloadButton"] button {{
            background-color: {t['btn_bg']} !important; color: {t['btn_text']} !important;
            border: 1px solid {t['btn_border']} !important; border-radius: 6px; font-weight: 500;
            transition: all 0.2s ease; box-shadow: 0 1px 2px rgba(0,0,0,0.05); text-decoration: none !important;
        }}
        .stButton button:hover, [data-testid="stLinkButton"] a:hover, [data-testid="stDownloadButton"] button:hover {{
            background-color: {t['btn_hover_bg']} !important; border-color: #A0AEC0 !important;
            transform: translateY(-1px); color: {t['text_main']} !important;
        }}

        /* 文本与排版 */
        .probe-title {{ font-size: 1.2rem; font-weight: 700; margin-bottom: 6px; color: {t['text_main']}; letter-spacing: -0.01em; }}
        .probe-meta {{ font-size: 0.9rem; color: {t['meta_color']}; font-family: 'Source Sans Pro', sans-serif; }}

        /* 徽章 */
        .badge-target {{ background-color: {t['badge_target_bg']}; color: {t['badge_target_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; font-weight: 600; display: inline-block; margin-right: 6px; }}
        .badge-type {{ background-color: {t['badge_type_bg']}; color: {t['badge_type_text']}; padding: 3px 10px; border-radius: 100px; font-size: 0.8rem; display: inline-block; }}
        .badge-new {{ background: linear-gradient(135deg, #FFD700 0%, #F59E0B 100%); color: white; padding: 2px 6px; border-radius: 4px; font-size: 0.7rem; font-weight: 800; margin-left: 8px; vertical-align: middle; box-shadow: 0 2px 5px rgba(245, 158, 11, 0.4); }}
    </style>
    """, unsafe_allow_html=True)

# ================= 4. 数据处理 =================
def load_data():
    if not os.path.exists(DATA_FILE):
        return pd.DataFrame()
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    return pd.DataFrame(data)

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

def extract_years(df):
    """自动提取最早和最晚年份"""
    if df.empty or 'date' not in df.columns:
        return "N/A", "N/A"
    try:
        # 使用正则提取4位数字年份
        years = df['date'].astype(str).str.extract(r'(\d{4})').astype(float)
        min_year = int(years.min().iloc[0])
        max_year = int(years.max().iloc[0])
        return min_year, max_year
    except:
        return "2024", datetime.datetime.now().year

# ================= 5. 渲染组件 =================

def render_sidebar(df, theme):
    with st.sidebar:
        # 1. 标题 (已修改)
        st.title("🧬 Auto Sensor Base")
        st.caption("Automated Tracking System")
        
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
        st.markdown("<br>", unsafe_allow_html=True)

        # 2. 下载按钮
        if not df.empty:
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Download CSV Dataset", csv, "probes_database.csv", "text/csv", use_container_width=True)
        
        st.divider()

        # 3. 筛选器
        filtered = df.copy()
        if not df.empty:
            df['target'] = df['target'].astype(str)
            df['color'] = df['color'].astype(str)
            
            sel_target = st.selectbox("Target Molecule", ["All"] + sorted(list(df['target'].unique())))
            sel_color = st.selectbox("Fluorescence Color", ["All"] + sorted(list(df['color'].unique())))
            
            if sel_target != "All": filtered = filtered[filtered['target'] == sel_target]
            if sel_color != "All": filtered = filtered[filtered['color'] == sel_color]

        # 4. 统计信息 (美化版)
        min_y, max_y = extract_years(filtered)
        st.markdown(f"""
        <div style='margin-top: 30px; padding: 15px; background: rgba(0,0,0,0.03); border-radius: 12px; text-align: center; border: 1px solid {theme['border_color']};'>
            <div style='font-size: 0.85rem; color: {theme['meta_color']}; text-transform: uppercase; letter-spacing: 1px;'>Total Sensors</div>
            <div style='font-size: 2rem; font-weight: 800; color: {theme['text_sidebar']}; line-height: 1.2;'>{len(filtered)}</div>
            <div style='font-size: 0.85rem; color: {theme['text_sidebar']}; opacity: 0.8; margin-top: 4px;'>
                ({min_y} - {max_y})
            </div>
        </div>
        """, unsafe_allow_html=True)

        # 5. 底部品牌信息 (Footer)
        st.markdown(f"""
        <div style='margin-top: 50px; padding-top: 20px; border-top: 1px solid {theme['border_color']}; text-align: center;'>
            <div style='font-weight: 600; font-size: 0.95rem; margin-bottom: 4px;'>Chimera Sensor Team</div>
            <div style='font-size: 0.8rem; opacity: 0.8; margin-bottom: 8px;'>Institute of Neuroscience, CAS</div>
            <div style='font-size: 0.85rem;'>
                <a href="http://www.cns.ac.cn" target="_blank" style='text-decoration: none; border-bottom: 1px dotted;'>www.cns.ac.cn</a>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        return is_light, filtered

def render_main_feed(df, theme):
    st.header("🚀 Latest Probes")

    if df.empty:
        st.info("No data available yet.")
        return

    try: df = df.sort_values(by='date', ascending=False)
    except: pass

    for index, row in df.iterrows():
        current_year = datetime.datetime.now().year
        pub_year = str(row.get('date', ''))
        is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
        new_badge = '<span class="badge-new">NEW</span>' if is_new else ""

        with st.container(border=True):
            c1, c2, c3 = st.columns([0.2, 5.5, 1.0])
            with c1: st.markdown(f"<div style='padding-top: 4px; text-align: center;'>{get_color_circle_html(row['color'])}</div>", unsafe_allow_html=True)
            with c2:
                st.markdown(f"""
                <div class="probe-title">{row['probe_name']} {new_badge}</div>
                <div style="margin-top: 8px; line-height: 1.8;">
                    <span class="badge-target">{row['target']}</span>
                    <span class="badge-type">{row.get('type', 'Unknown')}</span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{row.get('journal', 'Unknown')}</i></span>
                    <span style="color: {theme.get('border_color')}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {row.get('date', 'N/A')}</span>
                </div>
                """, unsafe_allow_html=True)
            with c3:
                st.markdown("<div style='height: 6px'></div>", unsafe_allow_html=True) 
                if row.get('doi') and "http" in row['doi']:
                    st.link_button("Read", row['doi'], use_container_width=True)
                else:
                    st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)

# ================= 6. 程序入口 =================
def main():
    df = load_data()
    # 先获取默认主题以渲染侧边栏
    is_light_mode = False 
    
    # 这里的逻辑稍微调整：因为 toggle 在侧边栏内部，我们需要先定义一个临时 theme 来渲染侧边栏框架
    # 但由于 Streamlit 的运行机制，我们可以直接在 render_sidebar 内部处理
    # 为了代码整洁，我们先用默认 Dark 渲染一次拿到 user input，然后全量渲染
    
    # 渲染侧边栏并获取用户的主题选择 + 筛选后的数据
    # 我们先传递一个临时的 Dark Theme 进去，因为此时还不知道用户选了啥，但侧边栏里的文字需要颜色
    temp_theme = get_theme_config(False) 
    
    # 注意：render_sidebar 需要 theme 参数来渲染底部的 Footer 颜色
    # 这是一个小小的 "先有鸡还是先有蛋" 问题。
    # 解决方法：我们把 toggle 放在最前面，根据 session state 或默认值拿到 is_light
    
    # --- 优化后的渲染顺序 ---
    # 1. 先加载数据
    # 2. 渲染 Sidebar 并直接获取 is_light 状态（Streamlit 会自动处理重运行）
    # 3. 根据状态生成 Config
    # 4. 注入 CSS
    # 5. 渲染 Main
    
    # 重新组织 render_sidebar 逻辑以支持动态传参
    with st.sidebar:
        st.title("🧬 Auto Sensor Base")
        st.caption("Automated Tracking System")
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
    
    theme_config = get_theme_config(is_light)
    inject_custom_css(theme_config)
    
    # 现在再次调用 sidebar 的剩余部分（筛选器、Footer等），传入正确的 theme
    _, filtered_df = render_sidebar_content(df, theme_config)
    
    render_main_feed(filtered_df, theme_config)

def render_sidebar_content(df, theme):
    # 这里只渲染侧边栏除 Title/Toggle 以外的内容
    with st.sidebar:
        st.markdown("<br>", unsafe_allow_html=True)
        if not df.empty:
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Download CSV Dataset", csv, "probes_database.csv", "text/csv", use_container_width=True)
        st.divider()

        filtered = df.copy()
        if not df.empty:
            df['target'] = df['target'].astype(str)
            df['color'] = df['color'].astype(str)
            sel_target = st.selectbox("Target Molecule", ["All"] + sorted(list(df['target'].unique())))
            sel_color = st.selectbox("Fluorescence Color", ["All"] + sorted(list(df['color'].unique())))
            if sel_target != "All": filtered = filtered[filtered['target'] == sel_target]
            if sel_color != "All": filtered = filtered[filtered['color'] == sel_color]

        min_y, max_y = extract_years(filtered)
        st.markdown(f"""
        <div style='margin-top: 30px; padding: 15px; background: rgba(0,0,0,0.03); border-radius: 12px; text-align: center; border: 1px solid {theme['border_color']};'>
            <div style='font-size: 0.85rem; color: {theme['meta_color']}; text-transform: uppercase; letter-spacing: 1px;'>Total Sensors</div>
            <div style='font-size: 2rem; font-weight: 800; color: {theme['text_sidebar']}; line-height: 1.2;'>{len(filtered)}</div>
            <div style='font-size: 0.85rem; color: {theme['text_sidebar']}; opacity: 0.8; margin-top: 4px;'>({min_y} - {max_y})</div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown(f"""
        <div style='margin-top: 50px; padding-top: 20px; border-top: 1px solid {theme['border_color']}; text-align: center;'>
            <div style='font-weight: 600; font-size: 0.95rem; margin-bottom: 4px; color: {theme['text_sidebar']};'>Chimera Nano Sensor Team</div>
            <div style='font-size: 0.8rem; opacity: 0.8; margin-bottom: 8px; color: {theme['text_sidebar']};'>Institute of Neuroscience, CAS</div>
            <div style='font-size: 0.85rem;'>
                <a href="http://www.cns.ac.cn" target="_blank" style='text-decoration: none; border-bottom: 1px dotted; color: {theme['text_sidebar']};'>www.cns.ac.cn</a>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        return None, filtered

if __name__ == "__main__":
    main()