import streamlit as st
import pandas as pd
from backend import COMPLETE_GENOMIC_DATABASE, process_uploaded_fasta, generate_coverage_plot

st.set_page_config(page_title="Lentera Gen - Decision Support", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: #f8fafc !important; color: #0f172a !important; font-family: 'Segoe UI', sans-serif; }
    label, p, span, h1, h2, h3, h4, .stMarkdown { color: #0f172a !important; }

    [data-testid="stSidebar"] {
        background-color: #ffffff !important;
        border-right: 2px solid #e2e8f0;
    }
    [data-testid="stSidebar"] label, [data-testid="stSidebar"] p {
        color: #0f172a !important;
        font-size: 15px !important;
        font-weight: 700 !important;
    }

    .header-banner {
        background: linear-gradient(135deg, #0a1b3d 0%, #2a52c9 100%);
        padding: 24px 30px;
        border-radius: 14px;
        margin-bottom: 24px;
        box-shadow: 0 4px 12px rgba(10, 27, 61, 0.15);
    }
    .header-title { font-size: 28px; font-weight: 800; color: #ffffff !important; margin-bottom: 6px; }
    .header-sub { font-size: 15px; color: #0abab5 !important; font-weight: 600; }
    
    .senior-metric-card {
        background-color: #ffffff !important;
        padding: 18px;
        border-radius: 12px;
        border-left: 6px solid #2a52c9;
        box-shadow: 0 2px 8px rgba(15, 23, 42, 0.06);
        margin-bottom: 16px;
    }
    .senior-metric-label { font-size: 13px; color: #475569 !important; font-weight: 700; text-transform: uppercase; }
    .senior-metric-value { font-size: 22px; font-weight: 800; color: #0a1b3d !important; margin-top: 4px; }

    .gene-detail-card {
        background-color: #ffffff !important;
        border-radius: 12px;
        padding: 20px;
        margin-bottom: 16px;
        border: 1px solid #cbd5e1;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
    }
    .badge-gene {
        background-color: #2a52c9 !important;
        color: #ffffff !important;
        padding: 6px 14px;
        border-radius: 8px;
        font-weight: 800;
        font-size: 16px;
    }
    .badge-efornas-green {
        background-color: #d1fae5 !important;
        color: #065f46 !important;
        padding: 6px 12px;
        border-radius: 8px;
        font-weight: 700;
        font-size: 13px;
        border: 1px solid #a7f3d0;
    }
    .badge-efornas-yellow {
        background-color: #fef3c7 !important;
        color: #92400e !important;
        padding: 6px 12px;
        border-radius: 8px;
        font-weight: 700;
        font-size: 13px;
        border: 1px solid #fde68a;
    }
    .badge-efornas-red {
        background-color: #ffe4e6 !important;
        color: #9f1239 !important;
        padding: 6px 12px;
        border-radius: 8px;
        font-weight: 700;
        font-size: 13px;
        border: 1px solid #fecdd3;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="header-banner">
    <div class="header-title">🧬 Lentera Gen — Clinical Decision Support Dashboard</div>
    <div class="header-sub">Panel Hotspot Kanker Payudara Kompak Terkalibrasi Genetik Indonesia & Verification E-FORNAS JKN</div>
</div>
""", unsafe_allow_html=True)

tab_panel, tab_coverage, tab_efornas = st.tabs([
    "🏆 Multi-Database Ranked Gene Panel", 
    "📈 Coverage Curve & Variant Calling",
    "💊 Panduan Pembiayaan E-FORNAS JKN"
])

# ==============================================================================
# TAB 1: RANKED PANEL (TERMASUK CLINICAL SIGNIFICANCE & ALLELE PAIR)
# ==============================================================================
with tab_panel:
    st.subheader("🏆 Multi-Database Ranked Target Gene Panel (13 Hotspots)")
    st.markdown("Pemeringkatan gen target terintegrasi dari CIViC, ClinVar, SMC Korea Cohort, Malaysian, dan Chinese Cohorts.")
    
    st.sidebar.header("🔍 Filter Panel Gen Target")
    selected_genes = st.sidebar.multiselect(
        "Pilih Gen Target:",
        options=["PIK3CA", "AKT1", "ERBB2", "PTEN", "BRCA2", "TP53"],
        default=["PIK3CA", "AKT1", "ERBB2", "PTEN", "BRCA2", "TP53"]
    )
    
    df_full = pd.DataFrame(COMPLETE_GENOMIC_DATABASE)
    df_filtered = df_full[df_full["Genes"].isin(selected_genes)]

    view_type = st.radio("Pilih Mode Tampilan:", ["Kartu Informasi Ringkas (Rekomendasi)", "Tabel Data Lengkap"], horizontal=True)
    
    if view_type == "Kartu Informasi Ringkas (Rekomendasi)":
        for idx, row in df_filtered.iterrows():
            efornas_str = row['Status_EFORNAS']
            if "🟢" in efornas_str:
                badge_style = "badge-efornas-green"
            elif "🟡" in efornas_str:
                badge_style = "badge-efornas-yellow"
            else:
                badge_style = "badge-efornas-red"
                
            st.markdown(f"""
            <div class="gene-detail-card">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div>
                        <span class="badge-gene">{row['Genes']} — {row['Hotspots']}</span>
                        <span style="margin-left: 12px; font-size: 16px; font-weight: 700; color: #2a52c9;">{row['HGVSp_Short']}</span>
                    </div>
                    <span class="{badge_style}">{row['Status_EFORNAS']}</span>
                </div>
                <div style="margin-top: 14px; font-size: 15px; color: #0f172a; line-height: 1.7;">
                    <strong>🎯 Rekomendasi Terapi Target:</strong> <span style="font-size: 16px; font-weight: 700; color: #0a1b3d;">{row['Matched Targeted Therapies']}</span><br>
                    <strong>🩺 Signifikansi Klinis:</strong> {row['Clinical Significance']} | <strong>Bukti Klinis:</strong> {row['Evidence Level']}<br>
                    <strong>🧬 Perubahan Alel (Ref ➔ Tumor):</strong> <code>{row['Reference_Allele']} ➔ {row['Tumor_Seq_Allele2']}</code> | <code>{row['HGVSc']}</code> (Chr {row['Chromosome']}:{row['Start_Position']})<br>
                    <small style="color: #64748b;">Sumber Database Alignment: {row['Source']}</small>
                </div>
            </div>
            """, unsafe_allow_html=True)
    else:
        st.dataframe(
            df_filtered[[
                "Genes", "Hotspots", "HGVSp_Short", "Clinical Significance",
                "Reference_Allele", "Tumor_Seq_Allele2", "Matched Targeted Therapies", 
                "Status_EFORNAS", "Evidence Level", "Source"
            ]],
            use_container_width=True,
            hide_index=True
        )

# ==============================================================================
# TAB 2: COVERAGE CURVE & ANALISIS VARIANT
# ==============================================================================
with tab_coverage:
    st.subheader("📈 Quality Control & Coverage Curve Pasien (9 Amplikon)")
    
    st.sidebar.markdown("---")
    st.sidebar.header("📁 Input Sekuens Pasien")
    
    input_mode = st.sidebar.radio("Pilih Sumber Sekuens:", ["Preset Sampel Demo", "Upload File FASTA Custom"])
    
    fasta_text = None
    if input_mode == "Preset Sampel Demo":
        sample_preset = st.sidebar.selectbox(
            "Pilih Sampel Pasien:",
            [
                "Sampel 1: PIK3CA H1047L Primary Hotspot",
                "Sampel 2: Low Coverage Alert (QC Fail Amp 4 & 8)",
                "Sampel 3: Triple Actionable Target (PIK3CA/ERBB2/TP53)",
                "Sampel 4: Kontrol Normal (Wild-Type)"
            ]
        )
        if "Sampel 1" in sample_preset: header_prefix = "SRR14149060"
        elif "Sampel 2" in sample_preset: header_prefix = "SRR14149061"
        elif "Sampel 3" in sample_preset: header_prefix = "SRR14149062"
        else: header_prefix = "SRR14149063"
            
        fasta_text = f">{header_prefix}.1 Amp_01_PIK3CA_H1047L length=144\nGCNTATGCACATTTTTAGAATTCCATTATCTTGAATACAAATAGGGAAATATCTTTTGGAAATTAAAATA\nATTTCAAAACAGTGGGGCAAAGAACCACTCACCTGTTGAATTTAACAGCAAGGCCCAGGTCAGCAATGCA\nGCAA\n"
    else:
        uploaded_file = st.sidebar.file_uploader("Upload File FASTA Pasien", type=["fasta", "fa", "txt"])
        if uploaded_file is not None:
            fasta_text = uploaded_file.getvalue().decode("utf-8")

    if fasta_text:
        results = process_uploaded_fasta(fasta_text)
        df_amp = results["amplicon_data"]
        
        st.info(f"📋 **Skenario Pasien Terdeteksi:** {results['scenario_title']}")
        
        mc1, mc2, mc3 = st.columns(3)
        with mc1:
            st.markdown(f'<div class="senior-metric-card"><div class="senior-metric-label">Total Read Sekuens</div><div class="senior-metric-value">{results["total_reads"]} Reads</div></div>', unsafe_allow_html=True)
        with mc2:
            st.markdown(f'<div class="senior-metric-card"><div class="senior-metric-label">Rata-rata Read Depth</div><div class="senior-metric-value" style="color:#2a52c9;">{results["mean_depth"]}x Coverage</div></div>', unsafe_allow_html=True)
        with mc3:
            qc_color = "#065f46" if "9/9" in results["passing_amplicons"] else "#e11d48"
            st.markdown(f'<div class="senior-metric-card"><div class="senior-metric-label">QC Amplikon Lolos (>=30x)</div><div class="senior-metric-value" style="color:{qc_color};">{results["passing_amplicons"]} PASS</div></div>', unsafe_allow_html=True)
            
        st.markdown("<br>", unsafe_allow_html=True)
        
        fig_plot = generate_coverage_plot(df_amp, results['scenario_title'])
        st.pyplot(fig_plot)
        
        st.markdown("<br>", unsafe_allow_html=True)
        st.subheader("📋 Status Alignment & Pemanggilan Mutasi per Amplikon")
        
        for idx, row in df_amp.iterrows():
            status_color = "#991b1b" if "MUTATED" in row["Hasil Alignment SNP"] else "#065f46"
            st.markdown(f"""
            <div style="background-color: #ffffff; padding: 14px 18px; border-radius: 10px; margin-bottom: 10px; border-left: 5px solid {status_color}; box-shadow: 0 1px 4px rgba(0,0,0,0.04);">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <strong>{row['Amplicon ID']} — Gen {row['Gen Target']} ({row['Hotspot SNP']})</strong>
                    <span style="color: {status_color}; font-weight: 800;">{row['Hasil Alignment SNP']}</span>
                </div>
                <div style="margin-top: 6px; font-size: 14px; color: #334155;">
                    Read Depth: <strong>{row['Read Depth (x)']}x</strong> | QC: <strong>{row['QC Threshold']}</strong> | Notasi: <code>{row['Notasi Genomik']}</code><br>
                    <strong>Therapy:</strong> {row['Rekomendasi Terapi']} | <strong>E-FORNAS:</strong> {row['Status E-FORNAS JKN']}
                </div>
            </div>
            """, unsafe_allow_html=True)
    else:
        st.info("👈 Silakan pilih preset sampel atau unggah file FASTA di sidebar untuk melihat analisis Coverage Curve.")

# ==============================================================================
# TAB 3: KATALOG PEMBIAYAAN E-FORNAS JKN
# ==============================================================================
with tab_efornas:
    st.subheader("💊 Panduan Pembiayaan Obat Target Kanker Payudara (E-FORNAS JKN)")
    st.markdown("Panduan ringkas bagi dokter klinis dan verifikator BPJS Kesehatan mengenai ketersediaan obat target di Formularium Nasional.")
    
    col_g, col_y, col_r = st.columns(3)
    
    with col_g:
        st.markdown("""
        <div style="background-color: #f0fdf4; border: 2px solid #bbf7d0; border-radius: 12px; padding: 18px;">
            <h4 style="color: #166534; margin: 0 0 10px 0;">🟢 Terdaftar E-FORNAS (Covered 100%)</h4>
            <ul style="font-size: 15px; color: #14532d; padding-left: 20px; line-height: 1.8;">
                <li><strong>Fulvestrant:</strong> Injeksi hormonal lini pertama/lanjutan.</li>
                <li><strong>Trastuzumab:</strong> Terapi anti-HER2 (Target mutasi ERBB2 V777L).</li>
                <li><strong>Doxorubicin:</strong> Kemoterapi standar (Target TP53).</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col_y:
        st.markdown("""
        <div style="background-color: #fffbeb; border: 2px solid #fde68a; border-radius: 12px; padding: 18px;">
            <h4 style="color: #92400e; margin: 0 0 10px 0;">🟡 Kombinasi Parsial JKN</h4>
            <ul style="font-size: 15px; color: #78350f; padding-left: 20px; line-height: 1.8;">
                <li><strong>Capivasertib + Fulvestrant:</strong> Fulvestrant ditanggung JKN (🟢), sedangkan Capivasertib belum terdaftar (🔴).</li>
                <li><strong>Alpelisib + Fulvestrant:</strong> Fulvestrant ditanggung JKN (🟢), Alpelisib belum terdaftar (🔴).</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col_r:
        st.markdown("""
        <div style="background-color: #fff1f2; border: 2px solid #fecdd3; border-radius: 12px; padding: 18px;">
            <h4 style="color: #9f1239; margin: 0 0 10px 0;">🔴 Belum Terdaftar E-FORNAS</h4>
            <ul style="font-size: 15px; color: #881337; padding-left: 20px; line-height: 1.8;">
                <li><strong>Alpelisib:</strong> PI3K Inhibitor.</li>
                <li><strong>Capivasertib:</strong> AKT Inhibitor.</li>
                <li><strong>Inavolisib:</strong> PI3Kalpha Inhibitor.</li>
                <li><strong>Neratinib:</strong> HER2 TKI Inhibitor.</li>
                <li><strong>Olaparib & Talazoparib:</strong> PARP Inhibitors (BRCA2).</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)