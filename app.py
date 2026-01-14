import streamlit as st
import pandas as pd
import json
import os

# ================= 页面配置 =================
st.set_page_config(
    page_title="🧬 遗传编码探针数据库",
    page_icon="🔬",
    layout="wide"
)

# 文件路径
DATA_FILE = "processed_probes.json"

# ================= 辅助函数 =================
def load_data():
    """读取 JSON 数据并转换为 DataFrame"""
    if not os.path.exists(DATA_FILE):
        return pd.DataFrame() # 返回空表
    
    with open(DATA_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    # 转换为 DataFrame 方便处理
    df = pd.DataFrame(data)
    return df

def color_badge(color_name):
    """根据荧光颜色返回不同的颜色点"""
    c = color_name.lower()
    if "green" in c: return "🟢"
    if "red" in c: return "🔴"
    if "blue" in c or "cyan" in c: return "🔵"
    if "yellow" in c: return "🟡"
    return "⚪"

# ================= 页面逻辑 =================

# 初始化 Session State (用于记录当前选了哪个探针)
if 'selected_probe_index' not in st.session_state:
    st.session_state.selected_probe_index = None

def go_back():
    """返回列表页"""
    st.session_state.selected_probe_index = None

# 1. 加载数据
df = load_data()

# 2. 侧边栏：标题与筛选
with st.sidebar:
    st.title("🔬 FP-Sensor Auto-DB")
    st.markdown("自动追踪最新的遗传编码荧光探针文献。")
    st.divider()
    
    if not df.empty:
        # 筛选器
        st.subheader("🔍 筛选")
        all_targets = ["全部"] + list(df['target'].unique())
        selected_target = st.selectbox("按检测底物筛选", all_targets)
        
        all_colors = ["全部"] + list(df['color'].unique())
        selected_color = st.selectbox("按颜色筛选", all_colors)
        
        # 应用筛选
        filtered_df = df.copy()
        if selected_target != "全部":
            filtered_df = filtered_df[filtered_df['target'] == selected_target]
        if selected_color != "全部":
            filtered_df = filtered_df[filtered_df['color'] == selected_color]
        
        st.info(f"共展示 {len(filtered_df)} 个探针")
    else:
        filtered_df = pd.DataFrame()
        st.warning("暂无数据，请先运行爬虫脚本。")

# ================= 主界面内容 =================

# 场景 A: 详情页 (如果用户点击了某个探针)
if st.session_state.selected_probe_index is not None:
    # 获取当前选中的行数据
    # 注意：这里需要从原始 df 获取，因为 index 是固定的
    try:
        probe = df.loc[st.session_state.selected_probe_index]
    except KeyError:
        st.session_state.selected_probe_index = None
        st.rerun()

    # ---- 详情页布局 ----
    st.button("← 返回列表", on_click=go_back)
    
    st.markdown(f"# {color_badge(probe['color'])} {probe['probe_name']}")
    st.caption(f"发表于 *{probe.get('journal', 'Unknown Journal')}* ({probe.get('date', 'Unknown Date')})")
    
    # 核心指标卡片
    col1, col2, col3, col4 = st.columns(4)
    with col1: st.metric("检测底物", probe['target'])
    with col2: st.metric("荧光颜色", probe['color'])
    with col3: st.metric("探针类型", probe.get('type', 'N/A'))
    with col4: 
        if probe.get('doi') and "http" in probe['doi']:
            st.link_button("🔗 阅读原文", probe['doi'])
        else:
            st.metric("DOI", "Unavailable")

    st.divider()
    
    st.subheader("📝 摘要")
    st.info(probe['abstract'])
    
    st.subheader("⚙️ 原始数据 (JSON)")
    st.json(probe.to_dict())

# 场景 B: 列表页 (默认展示)
else:
    st.title("🚀 最新发布的探针列表")
    
    if filtered_df.empty:
        st.info("👋 还没有找到新探针。请运行后台脚本抓取数据，或手动生成一些测试数据。")
    else:
        # 使用卡片式布局展示列表
        for index, row in filtered_df.iterrows():
            # 创建一个带边框的容器
            with st.container(border=True):
                c1, c2, c3 = st.columns([1, 4, 1])
                
                with c1:
                    # 显示大大的颜色图标
                    st.markdown(f"<h1 style='text-align: center;'>{color_badge(row['color'])}</h1>", unsafe_allow_html=True)
                
                with c2:
                    st.subheader(f"{row['probe_name']}")
                    st.markdown(f"**Target:** `{row['target']}` | **Type:** {row.get('type', 'Unknown')}")
                    st.markdown(f"*{row['title']}*")
                
                with c3:
                    st.markdown("<br>", unsafe_allow_html=True) # 占位符，为了按钮居中
                    # 点击按钮，更新 session_state，然后 rerun 刷新页面
                    if st.button("查看详情", key=f"btn_{index}"):
                        st.session_state.selected_probe_index = index
                        st.rerun()