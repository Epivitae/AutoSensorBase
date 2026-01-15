import streamlit as st
import pandas as pd
import json
import os
import datetime
import time
import requests
from github import Github
from Bio import Entrez

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
ZHIPU_API_KEY = st.secrets.get("ZHIPU_API_KEY") 
DATA_FILE = "processed_probes.json"

# 初始化 Biopython
Entrez.email = "wangk@ion.ac.cn" 

# ================= 2. 核心功能函数 =================

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
        repo.update_file(contents.path, "🤖 Admin: Manual data update", json_content, contents.sha)
        return True
    except Exception as e:
        st.error(f"GitHub Sync Error: {e}")
        return False

def fetch_pubmed_metadata(doi):
    """根据 DOI 从 PubMed 获取元数据"""
    try:
        # 1. 搜索 DOI 对应的 PMID
        search_handle = Entrez.esearch(db="pubmed", term=doi, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_record["IdList"]
        if not id_list:
            return None, "DOI not found in PubMed."

        # 2. 获取详细信息
        fetch_handle = Entrez.efetch(db="pubmed", id=id_list[0], rettype="medline", retmode="xml")
        data = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        if not data['PubmedArticle']:
             return None, "No article details found."

        article = data['PubmedArticle'][0]['MedlineCitation']['Article']
        
        # 提取字段
        title = article.get('ArticleTitle', 'No Title')
        journal = article.get('Journal', {}).get('Title', 'Unknown Journal')
        year = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', str(datetime.datetime.now().year))
        
        abstract_list = article.get('Abstract', {}).get('AbstractText', [])
        abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
        
        # === 🔧 修复点：自动补全 DOI 链接 ===
        final_doi = doi.strip()
        if not final_doi.startswith("http"):
            final_doi = f"https://doi.org/{final_doi}"

        return {
            "title": title,
            "journal": journal,
            "date": year,
            "abstract": abstract,
            "doi": final_doi # 现在是完整的 URL
        }, None
        
    except Exception as e:
        return None, str(e)

def ai_extract_attributes(abstract):
    """强制 AI 从摘要中提取属性"""
    if not ZHIPU_API_KEY:
        return {"target": "Manual", "color": "Other", "type": "Manual"}

    url = "https://open.bigmodel.cn/api/paas/v4/chat/completions"
    headers = {"Authorization": f"Bearer {ZHIPU_API_KEY}", "Content-Type": "application/json"}
    
    prompt = f"""
    Read this abstract and extract the fluorescent sensor's attributes.
    Abstract: {abstract}
    
    Return JSON ONLY with these keys:
    {{
        "target": "What molecule is detected? (e.g. Dopamine, Calcium)",
        "color": "Fluorescence color (Green, Red, Blue, Yellow, etc.)",
        "type": "Sensor type (e.g. CPFP, GPCR-based, snifit)"
    }}
    If unsure, use "Unknown". Do NOT return Markdown.
    """
    
    payload = {
        "model": "glm-4-flash",
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.1
    }

    try:
        resp = requests.post(url, headers=headers, json=payload, timeout=10)
        if resp.status_code == 200:
            content = resp.json()['choices'][0]['message']['content']
            clean_json = content.replace("```json", "").replace("```", "").strip()
            return json.loads(clean_json)
    except Exception as e:
        print(f"AI Error: {e}")
    
    return {"target": "Manual", "color": "Other", "type": "Manual"}

# ================= 3. 主题配色 =================
def get_theme_config(is_light_mode):
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
    st.markdown(f"""
    <style>
        .stApp {{ background-color: {t['bg_color']}; color: {t['text_main']}; }}
        .stDecoration {{ display: none !important; }}
        header[data-testid="stHeader"] {{ background-color: transparent !important; visibility: {t['header_visibility']}; }}
        [data-testid="stSidebar"] {{ background-color: {t['sidebar_bg']}; border-right: 1px solid #E2E8F0; }}
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] span, 
        [data-testid="stSidebar"] div, [data-testid="stSidebar"] label, [data-testid="stSidebar"] p {{ color: {t['text_sidebar']} !important; }}
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

# ================= 4. 数据处理 =================
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

# ================= 5. 管理员面板 (Smart Add) =================
def render_admin_panel(current_data):
    with st.sidebar:
        st.markdown("---")
        st.markdown("### 🛠️ Admin Panel (AI Auto-Fill)")
        
        with st.expander("✨ Smart Add Probe", expanded=True):
            with st.form("smart_add_form"):
                input_doi = st.text_input("1. DOI (Required)", placeholder="e.g. 10.1038/s41592-024-02202-x")
                input_name = st.text_input("2. Probe Name (Required)", placeholder="e.g. GCaMP8")
                
                st.caption("Leave the rest to AI! 🤖")
                submitted = st.form_submit_button("🚀 Fetch & Add")
                
                if submitted:
                    if not input_doi or not input_name:
                        st.error("Please fill in DOI and Name.")
                    else:
                        status_msg = st.empty()
                        status_msg.info("⏳ Fetching from PubMed...")
                        
                        # 1. 爬取 PubMed
                        meta_data, error = fetch_pubmed_metadata(input_doi)
                        
                        if error:
                            status_msg.error(f"PubMed Error: {error}")
                        else:
                            status_msg.info("🧠 AI Extracting Attributes...")
                            
                            # 2. AI 提取属性
                            ai_attrs = ai_extract_attributes(meta_data['abstract'])
                            
                            # 3. 合并数据
                            new_entry = {
                                "probe_name": input_name,
                                "is_new": True,
                                **meta_data,
                                **ai_attrs
                            }
                            
                            # 4. 保存
                            current_data.insert(0, new_entry)
                            if update_github_data(current_data):
                                status_msg.success(f"✅ Added {input_name} successfully!")
                                time.sleep(1)
                                st.rerun()

# ================= 6. 渲染侧边栏 =================
def render_sidebar_content(data, theme):
    with st.sidebar:
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

        min_y, max_y = extract_years(filtered)
        st.markdown(f"""
        <div style='margin-top: 20px; padding: 15px; background: rgba(0,0,0,0.03); border-radius: 12px; text-align: center; border: 1px solid {theme['border_color']};'>
            <div style='font-size: 0.8rem; color: {theme['meta_color']}; letter-spacing: 1px;'>TOTAL SENSORS</div>
            <div style='font-size: 2rem; font-weight: 800; color: {theme['text_sidebar']};'>{len(filtered)}</div>
            <div style='font-size: 0.8rem; color: {theme['text_sidebar']}; opacity: 0.8;'>({min_y} - {max_y})</div>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("<br>", unsafe_allow_html=True)
        
        # 登录逻辑
        if 'is_admin' not in st.session_state: st.session_state.is_admin = False
        
        with st.expander("🔐 Admin Login", expanded=st.session_state.is_admin):
            if not st.session_state.is_admin:
                with st.form("login_form"):
                    pwd = st.text_input("Password", type="password")
                    if st.form_submit_button("Login"):
                        if pwd == ADMIN_PWD:
                            st.session_state.is_admin = True
                            st.rerun()
                        else:
                            st.error("Wrong password")
            else:
                st.success("✅ Logged in")
                if st.button("Logout"):
                    st.session_state.is_admin = False
                    st.rerun()
        
        st.markdown(f"""
        <div style='margin-top: 40px; padding-top: 20px; border-top: 1px solid {theme['border_color']}; text-align: center;'>
            <div style='font-weight: 600; font-size: 0.9rem; margin-bottom: 4px; color: {theme['text_sidebar']};'>Chimera Nano Sensor Team</div>
            <div style='font-size: 0.8rem; opacity: 0.8; margin-bottom: 8px; color: {theme['text_sidebar']};'>Institute of Neuroscience, CAS</div>
            <div style='font-size: 0.8rem;'>
                <a href="http://www.cns.ac.cn" target="_blank" style='text-decoration: none; border-bottom: 1px dotted; color: {theme['text_sidebar']};'>www.cns.ac.cn</a>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        return filtered, st.session_state.is_admin

# ================= 7. 渲染主列表 =================
def render_main_feed(data, theme, is_admin):
    st.header("Latest Genetically Encoded Fluorescent Sensors")
    if not data:
        st.info("No data available.")
        return

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
                if is_admin:
                    if st.button("🗑️ Delete", key=f"del_{index}", type="primary", use_container_width=True):
                        data.pop(index)
                        if update_github_data(data):
                            st.success("Deleted!"); import time; time.sleep(0.5); st.rerun()
                else:
                    if row.get('doi') and "http" in row['doi']:
                        st.link_button("Read", row['doi'], use_container_width=True)
                    else:
                        st.button("No Link", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            with st.expander("View Abstract", expanded=False):
                st.markdown(f"<div style='opacity: 0.85; line-height: 1.6;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)

# ================= 8. 程序入口 =================
def main():
    data_list = load_data()
    
    with st.sidebar:
        st.title("Auto Sensor Base")
        st.caption("Automated Tracking System")
        is_light = st.toggle("🌞 Light Mode / 🌜 Dark", value=False, key="theme_toggle")

    theme_config = get_theme_config(is_light)
    inject_custom_css(theme_config)
    
    filtered_data, is_admin = render_sidebar_content(data_list, theme_config)
    
    # 如果是管理员，显示智能添加面板
    if is_admin:
        render_admin_panel(data_list)
    
    render_main_feed(filtered_data, theme_config, is_admin)

if __name__ == "__main__":
    main()