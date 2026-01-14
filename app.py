import streamlit as st
import pandas as pd
import json
import os
import datetime

# ================= 1. 配置与初始化 =================
st.set_page_config(
    page_title="FP-Sensor Auto-DB",
    page_icon="🧬",
    layout="wide"
)

DATA_FILE = "processed_probes.json"

# ================= 2. 主题配色定义 (方便后期修改) =================
def get_theme_config(is_light_mode):
    """
    根据模式返回配色字典。
    以后如果想改颜色，只需要改这里的代码即可。
    """
    if is_light_mode:
        # --- 🌞 白天模式 (Light Mode) ---
        return {
            "mode": "light",
            "bg_color": "#F7F9FB",          # 全局背景：极淡灰蓝
            "sidebar_bg": "#FFFFFF",        # 侧边栏：纯白
            "text_main": "#2D3748",         # 主文字：深灰
            "text_sidebar": "#1A202C",      # 侧边栏文字：深黑
            "meta_color": "#718096",        # 辅助文字：中灰
            
            # 卡片
            "card_bg": "#FFFFFF",
            "card_border": "1px solid #E2E8F0",
            "card_shadow": "0 2px 5px rgba(0,0,0,0.03)",
            
            # 按钮 (通用)
            "btn_bg": "#FFFFFF",
            "btn_text": "#4A5568",
            "btn_border": "#CBD5E0",
            "btn_hover_bg": "#F7FAFC",
            
            # 徽章 (Tag)
            "badge_target_bg": "#E6FFFA",   # 薄荷绿背景
            "badge_target_text": "#2C7A7B", # 深青色文字
            "badge_type_bg": "#EDF2F7",     # 浅灰背景
            "badge_type_text": "#4A5568",   # 深灰文字
            
            "header_visibility": "hidden"   # 隐藏默认黑头
        }
    else:
        # --- 🌜 黑夜模式 (Dark Mode) ---
        return {
            "mode": "dark",
            "bg_color": "#0E1117",
            "sidebar_bg": "#262730",
            "text_main": "#FAFAFA",
            "text_sidebar": "#FAFAFA",
            "meta_color": "#9CA3AF",
            
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
    """
    将 Python 字典转换为 CSS 样式并注入页面
    """
    st.markdown(f"""
    <style>
        /* [全局] 背景与文字 */
        .stApp {{
            background-color: {t['bg_color']};
            color: {t['text_main']};
        }}
        
        /* [顶部] 隐藏装饰条与 Header */
        .stDecoration {{ display: none !important; }}
        header[data-testid="stHeader"] {{
            background-color: transparent !important;
            visibility: {t['header_visibility']};
        }}
        
        /* [侧边栏] 强制样式覆盖 */
        [data-testid="stSidebar"] {{
            background-color: {t['sidebar_bg']};
            border-right: 1px solid #E2E8F0;
        }}
        /* 强制侧边栏所有元素文字颜色 (解决白天模式看不清的问题) */
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3,
        [data-testid="stSidebar"] span, [data-testid="stSidebar"] label, [data-testid="stSidebar"] div {{
            color: {t['text_sidebar']} !important;
        }}
        /* 侧边栏输入框背景修复 */
        [data-testid="stSidebar"] .stSelectbox > div > div {{
            background-color: {t['bg_color']};
            color: {t['text_main']};
            border-color: {t['btn_border']};
        }}

        /* [布局] 调整顶部间距 */
        .block-container {{ padding-top: 1.5rem; padding-bottom: 3rem; }}

        /* [卡片] 样式 */
        [data-testid="stVerticalBlockBorderWrapper"] > div {{
            background-color: {t['card_bg']};
            border: {t['card_border']} !important;
            box-shadow: {t['card_shadow']};
            border-radius: 10px;
            padding: 1.2rem;
            transition: box-shadow 0.2s ease;
        }}
        
        /* [按钮] 万能覆盖 (包括 Download, LinkButton, Button) */
        .stButton button, 
        [data-testid="stLinkButton"] a, 
        [data-testid="stDownloadButton"] button {{
            background-color: {t['btn_bg']} !important;
            color: {t['btn_text']} !important;
            border: 1px solid {t['btn_border']} !important;
            border-radius: 6px;
            font-weight: 500;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(0,0,0,0.05);
            text-decoration: none !important;
        }}
        
        /* 按钮悬停态 */
        .stButton button:hover, 
        [data-testid="stLinkButton"] a:hover, 
        [data-testid="stDownloadButton"] button:hover {{
            background-color: {t['btn_hover_bg']} !important;
            border-color: #A0AEC0 !important;
            transform: translateY(-1px);
            color: {t['text_main']} !important;
        }}

        /* [排版] 标题与元数据 */
        .probe-title {{
            font-size: 1.2rem; font-weight: 700; margin-bottom: 6px;
            color: {t['text_main']}; letter-spacing: -0.01em;
        }}
        .probe-meta {{
            font-size: 0.9rem; color: {t['meta_color']};
            font-family: 'Source Sans Pro', sans-serif;
        }}

        /* [徽章] 胶囊样式 */
        .badge-target {{
            background-color: {t['badge_target_bg']};
            color: {t['badge_target_text']};
            padding: 3px 10px; border-radius: 100px;
            font-size: 0.8rem; font-weight: 600;
            display: inline-block; margin-right: 6px;
        }}
        .badge-type {{
            background-color: {t['badge_type_bg']};
            color: {t['badge_type_text']};
            padding: 3px 10px; border-radius: 100px;
            font-size: 0.8rem; display: inline-block;
        }}
        
        /* NEW 星标 */
        .badge-new {{
            background: linear-gradient(135deg, #FFD700 0%, #F59E0B 100%);
            color: white; padding: 2px 6px; border-radius: 4px;
            font-size: 0.7rem; font-weight: 800; margin-left: 8px;
            vertical-align: middle; box-shadow: 0 2px 5px rgba(245, 158, 11, 0.4);
        }}
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
    """生成带颜色的小圆点 HTML"""
    c = str(color_name).lower()
    hex_color = "#D1D5DB"
    if "green" in c: hex_color = "#10B981"
    elif "red" in c: hex_color = "#EF4444"
    elif "blue" in c or "cyan" in c: hex_color = "#3B82F6"
    elif "yellow" in c or "gold" in c: hex_color = "#F59E0B"
    elif "orange" in c: hex_color = "#F97316"
    elif "purple" in c: hex_color = "#8B5CF6"
    
    return f"""
    <div style="width: 14px; height: 14px; background-color: {hex_color}; 
    border-radius: 50%; display: inline-block; box-shadow: 0 0 0 2px {hex_color}30;"></div>
    """

# ================= 5. 渲染组件 =================

def render_sidebar(df):
    """渲染侧边栏"""
    with st.sidebar:
        st.title("🧬 Auto-DB")
        st.caption("Automated Tracking System")
        
        # 主题切换
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False)
        
        st.markdown("<br>", unsafe_allow_html=True)

        # 下载按钮
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

        # 筛选器
        selected_target = "All"
        selected_color = "All"
        
        if not df.empty:
            df['target'] = df['target'].astype(str)
            df['color'] = df['color'].astype(str)
            
            all_targets = ["All"] + sorted(list(df['target'].unique()))
            selected_target = st.selectbox("Target Molecule", all_targets)
            
            all_colors = ["All"] + sorted(list(df['color'].unique()))
            selected_color = st.selectbox("Fluorescence Color", all_colors)
            
            # 底部统计
            filtered = df.copy()
            if selected_target != "All": filtered = filtered[filtered['target'] == selected_target]
            if selected_color != "All": filtered = filtered[filtered['color'] == selected_color]
            
            st.markdown(f"""
            <div style='margin-top: 20px; text-align: center; opacity: 0.7;'>
                Found <b>{len(filtered)}</b> probes
            </div>
            """, unsafe_allow_html=True)
            
            return is_light, filtered
        
        return is_light, pd.DataFrame()

def render_main_feed(df, theme):
    """渲染主列表"""
    st.header("🚀 Latest Probes")

    if df.empty:
        st.info("No data available yet.")
        return

    # 按日期排序
    try: df = df.sort_values(by='date', ascending=False)
    except: pass

    for index, row in df.iterrows():
        # 新探针标记逻辑 (当年或明年)
        current_year = datetime.datetime.now().year
        pub_year = str(row.get('date', ''))
        is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
        new_badge = '<span class="badge-new">NEW</span>' if is_new else ""

        # --- 卡片容器 ---
        with st.container(border=True):
            c1, c2, c3 = st.columns([0.2, 5.5, 1.0])
            
            with c1: # 颜色点
                st.markdown(f"<div style='padding-top: 4px; text-align: center;'>{get_color_circle_html(row['color'])}</div>", unsafe_allow_html=True)
            
            with c2: # 主要信息
                st.markdown(f"""
                <div class="probe-title">{row['probe_name']} {new_badge}</div>
                <div style="margin-top: 8px; line-height: 1.8;">
                    <span class="badge-target">{row['target']}</span>
                    <span class="badge-type">{row.get('type', 'Unknown')}</span>
                    <span style="color: {theme.get('btn_border')}; margin: 0 8px;">|</span>
                    <span class="probe-meta"><i>{row.get('journal', 'Unknown')}</i></span>
                    <span style="color: {theme.get('btn_border')}; margin: 0 8px;">•</span>
                    <span class="probe-meta">📅 {row.get('date', 'N/A')}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3: # 按钮
                st.markdown("<div style='height: 6px'></div>", unsafe_allow_html=True) 
                if row.get('doi') and "http" in row['doi']:
                    st.link_button("Read", row['doi'], use_container_width=True)
                else:
                    st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # 摘要
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)

# ================= 6. 程序入口 =================
def main():
    # 1. 加载数据
    df = load_data()
    
    # 2. 渲染侧边栏并获取状态
    is_light_mode, filtered_df = render_sidebar(df)
    
    # 3. 获取主题配置
    theme_config = get_theme_config(is_light_mode)
    
    # 4. 注入样式
    inject_custom_css(theme_config)
    
    # 5. 渲染主内容
    render_main_feed(filtered_df, theme_config)

if __name__ == "__main__":
    main()