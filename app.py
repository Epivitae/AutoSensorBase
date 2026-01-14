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

# Custom CSS for compact layout and badges
st.markdown("""
<style>
    /* Reduce vertical padding in block containers */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }
    /* Compact text styles */
    .probe-title {
        font-size: 1.2rem;
        font-weight: bold;
        margin-bottom: 0px;
    }
    .probe-meta {
        font-size: 0.85rem;
        color: #666;
        margin-bottom: 4px;
    }
    /* Badges */
    .badge-target {
        background-color: #e8fdf5;
        color: #0c8558;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.75rem;
        font-weight: 600;
        border: 1px solid #b7ebd6;
    }
    .badge-type {
        background-color: #f0f2f6;
        color: #31333F;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.75rem;
        border: 1px solid #dce0e6;
    }
    /* "New" Star Badge */
    .badge-new {
        background: linear-gradient(45deg, #FFD700, #FFC107);
        color: #7a5c00;
        padding: 2px 6px;
        border-radius: 4px;
        font-size: 0.7rem;
        font-weight: bold;
        margin-left: 8px;
        box-shadow: 0 1px 2px rgba(0,0,0,0.1);
        display: inline-block;
        vertical-align: middle;
    }
    /* Button compacting */
    .stButton button {
        height: 2.0rem;
        padding-top: 0;
        padding-bottom: 0;
    }
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
    """Returns a colored SVG circle for better visuals than emojis"""
    c = str(color_name).lower()
    hex_color = "#ccc" # default gray
    if "green" in c: hex_color = "#2ecc71"
    elif "red" in c: hex_color = "#e74c3c"
    elif "blue" in c or "cyan" in c: hex_color = "#3498db"
    elif "yellow" in c or "gold" in c: hex_color = "#f1c40f"
    elif "orange" in c: hex_color = "#e67e22"
    elif "purple" in c: hex_color = "#9b59b6"
    
    return f"""
    <div style="
        width: 18px; 
        height: 18px; 
        background-color: {hex_color}; 
        border-radius: 50%; 
        display: inline-block;
        border: 2px solid rgba(0,0,0,0.1);
        vertical-align: middle;
    "></div>
    """

# Load Data
df = load_data()

# ================= Sidebar =================
with st.sidebar:
    st.title("🧬 Auto-DB")
    st.caption("Genetically Encoded Fluorescent Probes")
    
    # Download Button
    if not df.empty:
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download CSV",
            data=csv,
            file_name='probes_database.csv',
            mime='text/csv',
        )
    
    st.divider()

    # Filters
    if not df.empty:
        # Normalize data types
        df['target'] = df['target'].astype(str)
        df['color'] = df['color'].astype(str)
        df['date'] = df['date'].astype(str)
        
        # Sort filters
        all_targets = ["All"] + sorted(list(df['target'].unique()))
        selected_target = st.selectbox("Filter by Target", all_targets)
        
        all_colors = ["All"] + sorted(list(df['color'].unique()))
        selected_color = st.selectbox("Filter by Color", all_colors)
        
        # Apply Filters
        filtered_df = df.copy()
        if selected_target != "All":
            filtered_df = filtered_df[filtered_df['target'] == selected_target]
        if selected_color != "All":
            filtered_df = filtered_df[filtered_df['color'] == selected_color]
        
        st.markdown(f"**Showing {len(filtered_df)} probes**")
    else:
        filtered_df = pd.DataFrame()

# ================= Main Content =================
st.header("🚀 Latest Probes")

if filtered_df.empty:
    st.info("No data available yet. Please run the backend script.")
else:
    # Sort by date (newest first) assuming date is roughly parsable
    # If date is just Year (2025), this works. 
    try:
        filtered_df = filtered_df.sort_values(by='date', ascending=False)
    except:
        pass

    # Loop through rows
    for index, row in filtered_df.iterrows():
        # Determine if "New" (Current year or next year)
        # Adjust '2025' to whatever logic you prefer
        current_year = datetime.datetime.now().year
        pub_year = str(row.get('date', ''))
        
        # Mark as new if published in current year or next year (preprint)
        is_new = str(current_year) in pub_year or str(current_year + 1) in pub_year
        new_badge = '<span class="badge-new">⭐ NEW</span>' if is_new else ""

        # Layout: Compact Container
        with st.container(border=True):
            # Columns: [Color Indicator] [Main Info] [Action Button]
            c1, c2, c3 = st.columns([0.3, 5, 1])
            
            with c1:
                # Vertical align the color dot
                st.markdown(f"<div style='margin-top: 5px;'>{color_circle(row['color'])}</div>", unsafe_allow_html=True)
            
            with c2:
                # Title + New Badge
                st.markdown(f"""
                <div class="probe-title">
                    {row['probe_name']} {new_badge}
                </div>
                """, unsafe_allow_html=True)
                
                # Meta info (Target | Type | Journal | Date) in one compact line
                target = row['target']
                ptype = row.get('type', 'Unknown')
                journal = row.get('journal', 'Unknown Journal')
                date = row.get('date', 'N/A')
                
                st.markdown(f"""
                <div style="margin-top: 4px;">
                    <span class="badge-target">{target}</span>
                    <span class="badge-type">{ptype}</span>
                    <span style="color: #bbb; margin: 0 6px;">|</span>
                    <span class="probe-meta"><i>{journal}</i></span>
                    <span style="color: #bbb; margin: 0 6px;">•</span>
                    <span class="probe-meta">📅 {date}</span>
                </div>
                """, unsafe_allow_html=True)

            with c3:
                # DOI Link Button
                if row.get('doi') and "http" in row['doi']:
                    st.link_button("Read", row['doi'], use_container_width=True)
                else:
                    st.button("No DOI", disabled=True, key=f"btn_{index}", use_container_width=True)
            
            # Abstract Expander (Compact)
            with st.expander(f"View Abstract: {row['title']}", expanded=False):
                st.markdown(f"<div style='font-size: 0.9rem; color: #444;'>{row.get('abstract', 'No abstract')}</div>", unsafe_allow_html=True)