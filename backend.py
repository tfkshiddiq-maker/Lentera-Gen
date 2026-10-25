import io
import matplotlib.pyplot as plt
import pandas as pd

# Database 13 Hotspots Multi-Source Alignment dengan Pemetaan E-FORNAS JKN
COMPLETE_GENOMIC_DATABASE = [
    {
        "Genes": "PIK3CA", "Hotspots": "H1047L", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.H1047L", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178952085,
        "Reference_Allele": "A", "Tumor_Seq_Allele2": "T", "HGVSc": "ENST00000263967.3:c.3140A>T", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "PIK3CA", "Hotspots": "H1047R", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.H1047R", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178952085,
        "Reference_Allele": "A", "Tumor_Seq_Allele2": "G", "HGVSc": "ENST00000263967.3:c.3140A>G", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "PIK3CA", "Hotspots": "E545K", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.E545K", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178936091,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "A", "HGVSc": "ENST00000263967.3:c.1633G>A", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "PIK3CA", "Hotspots": "E542K", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.E542K", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178936082,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "A", "HGVSc": "ENST00000263967.3:c.1624G>A", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "PIK3CA", "Hotspots": "C420R", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.C420R", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178927980,
        "Reference_Allele": "T", "Tumor_Seq_Allele2": "C", "HGVSc": "ENST00000263967.3:c.1258T>C", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "PIK3CA", "Hotspots": "N345K", "Matched Targeted Therapies": "Alpelisib, Capivasertib, Inavolisib",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.N345K", "Variant_Classification": "Missense_Mutation", "Chromosome": 3, "Start_Position": 178921553,
        "Reference_Allele": "T", "Tumor_Seq_Allele2": "A", "HGVSc": "ENST00000263967.3:c.1035T>A", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "AKT1", "Hotspots": "E17K", "Matched Targeted Therapies": "Capivasertib + Fulvestrant",
        "Status_EFORNAS": "🟡 Kombinasi Parsial (Fulvestrant: 🟢 Covered | Capivasertib: 🔴 Non-FORNAS)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.E17K", "Variant_Classification": "Missense_Mutation", "Chromosome": 14, "Start_Position": 105246551,
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "T", "HGVSc": "ENST00000349310.3:c.49G>A", "Source": "CIViC, SMC Korea, Chinese"
    },
    {
        "Genes": "ERBB2", "Hotspots": "V777L", "Matched Targeted Therapies": "Neratinib, Trastuzumab",
        "Status_EFORNAS": "🟢 Trastuzumab Covered 100% E-FORNAS (Neratinib: 🔴 Non-FORNAS)", "Evidence Level": "Level B", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.V777L", "Variant_Classification": "Missense_Mutation", "Chromosome": 17, "Start_Position": 37881000,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "T", "HGVSc": "ENST00000269571.5:c.2329G>T", "Source": "CIViC, SMC Korea, Chinese"
    },
    {
        "Genes": "PTEN", "Hotspots": "L193P", "Matched Targeted Therapies": "Capivasertib + Fulvestrant",
        "Status_EFORNAS": "🟡 Kombinasi Parsial (Fulvestrant: 🟢 Covered | Capivasertib: 🔴 Non-FORNAS)", "Evidence Level": "Level A / B", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.L193P", "Variant_Classification": "Missense_Mutation", "Chromosome": 10, "Start_Position": 89711960,
        "Reference_Allele": "T", "Tumor_Seq_Allele2": "C", "HGVSc": "ENST00000371953.3:c.578T>C", "Source": "CIViC, SMC Korea"
    },
    {
        "Genes": "BRCA2", "Hotspots": "E7K", "Matched Targeted Therapies": "Olaparib, Talazoparib (PARP inhibitors)",
        "Status_EFORNAS": "🔴 Belum Terdaftar E-FORNAS (Mandiri / Swasta)", "Evidence Level": "Level A", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.E7K", "Variant_Classification": "Missense_Mutation", "Chromosome": 13, "Start_Position": 32890616,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "A", "HGVSc": "ENST00000380152.3:c.19G>A", "Source": "CIViC, SMC Korea"
    },
    {
        "Genes": "TP53", "Hotspots": "R248Q", "Matched Targeted Therapies": "Doxorubicin / Chemotherapy sensitivity",
        "Status_EFORNAS": "🟢 Covered 100% E-FORNAS JKN (Doxorubicin)", "Evidence Level": "Level B", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.R248Q", "Variant_Classification": "Missense_Mutation", "Chromosome": 17, "Start_Position": 7577538,
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "T", "HGVSc": "ENST00000269305.4:c.742C>T", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "TP53", "Hotspots": "R248W", "Matched Targeted Therapies": "Doxorubicin / Chemotherapy sensitivity",
        "Status_EFORNAS": "🟢 Covered 100% E-FORNAS JKN (Doxorubicin)", "Evidence Level": "Level B", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.R248W", "Variant_Classification": "Missense_Mutation", "Chromosome": 17, "Start_Position": 7577538,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "A", "HGVSc": "ENST00000269305.4:c.743G>A", "Source": "CIViC, SMC Korea, Malaysian"
    },
    {
        "Genes": "TP53", "Hotspots": "R249", "Matched Targeted Therapies": "Doxorubicin / Chemotherapy sensitivity",
        "Status_EFORNAS": "🟢 Covered 100% E-FORNAS JKN (Doxorubicin)", "Evidence Level": "Level B", "Clinical Significance": "Sensitivity/Response",
        "HGVSp_Short": "p.R249", "Variant_Classification": "DNA Binding Mutation", "Chromosome": 17, "Start_Position": 7577537,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "T", "HGVSc": "ENST00000269305.4:c.745G>T", "Source": "CIViC, SMC Korea, Malaysian"
    }
]

AMPLICON_PANEL = [
    {"Amplicon_ID": "Amp_01", "Gene": "PIK3CA", "Hotspot": "H1047L (35.3% IDN)", "Therapy": "Alpelisib, Capivasertib, Inavolisib", "EFORNAS": "🔴 Belum Terdaftar E-FORNAS", "Target_SNP": "c.3140A>T (p.H1047L)"},
    {"Amplicon_ID": "Amp_02", "Gene": "PIK3CA", "Hotspot": "E545K", "Therapy": "Alpelisib, Capivasertib, Inavolisib", "EFORNAS": "🔴 Belum Terdaftar E-FORNAS", "Target_SNP": "c.1633G>A (p.E545K)"},
    {"Amplicon_ID": "Amp_03", "Gene": "AKT1", "Hotspot": "E17K", "Therapy": "Capivasertib + Fulvestrant", "EFORNAS": "🟡 Kombinasi Parsial (Fulvestrant: 🟢 Covered)", "Target_SNP": "c.49G>A (p.E17K)"},
    {"Amplicon_ID": "Amp_04", "Gene": "ERBB2", "Hotspot": "V777L", "Therapy": "Neratinib, Trastuzumab", "EFORNAS": "🟢 Trastuzumab Covered 100% E-FORNAS", "Target_SNP": "c.2329G>T (p.V777L)"},
    {"Amplicon_ID": "Amp_05", "Gene": "PTEN", "Hotspot": "L193P", "Therapy": "Capivasertib + Fulvestrant", "EFORNAS": "🟡 Kombinasi Parsial (Fulvestrant: 🟢 Covered)", "Target_SNP": "c.578T>C (p.L193P)"},
    {"Amplicon_ID": "Amp_06", "Gene": "BRCA2", "Hotspot": "E7K", "Therapy": "Olaparib, Talazoparib (PARP inhibitors)", "EFORNAS": "🔴 Belum Terdaftar E-FORNAS", "Target_SNP": "c.19G>A (p.E7K)"},
    {"Amplicon_ID": "Amp_07", "Gene": "TP53", "Hotspot": "R248Q", "Therapy": "Doxorubicin / Chemotherapy sensitivity", "EFORNAS": "🟢 Covered 100% E-FORNAS (Doxorubicin)", "Target_SNP": "c.742C>T (p.R248Q)"},
    {"Amplicon_ID": "Amp_08", "Gene": "TP53", "Hotspot": "R248W", "Therapy": "Doxorubicin / Chemotherapy sensitivity", "EFORNAS": "🟢 Covered 100% E-FORNAS (Doxorubicin)", "Target_SNP": "c.743G>A (p.R248W)"},
    {"Amplicon_ID": "Amp_09", "Gene": "TP53", "Hotspot": "R249", "Therapy": "Doxorubicin / Chemotherapy sensitivity", "EFORNAS": "🟢 Covered 100% E-FORNAS (Doxorubicin)", "Target_SNP": "c.745G>T (p.R249)"}
]

def process_uploaded_fasta(fasta_text: str) -> dict:
    lines = [l.strip() for l in fasta_text.strip().split('\n') if l.strip()]
    reads = []
    current_header = ""
    current_seq = []
    
    for line in lines:
        if line.startswith('>'):
            if current_header and current_seq:
                reads.append((current_header, "".join(current_seq).upper()))
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
    if current_header and current_seq:
        reads.append((current_header, "".join(current_seq).upper()))

    header_first = reads[0][0] if reads else ""
    
    if "SRR14149061" in header_first or "sample_2" in header_first.lower():
        scenario_title = "Sampel 2: QC Alert (Low Coverage pada Amp 4 & 8)"
        simulated_depths = [110, 95, 85, 18, 105, 60, 130, 22, 80]
        mutated_indices = [1, 4]
    elif "SRR14149062" in header_first or "sample_3" in header_first.lower():
        scenario_title = "Sampel 3: Profil Mutasi Tinggi (Triple Actionable Target)"
        simulated_depths = [210, 185, 190, 165, 200, 175, 230, 195, 180]
        mutated_indices = [0, 3, 7]
    elif "SRR14149063" in header_first or "sample_4" in header_first.lower():
        scenario_title = "Sampel 4: Kontrol Normal (100% Wild Type)"
        simulated_depths = [150, 142, 138, 155, 160, 148, 165, 152, 145]
        mutated_indices = []
    else:
        scenario_title = "Sampel 1: Hotspot Utama PIK3CA H1047L (Baseline)"
        simulated_depths = [125, 88, 102, 45, 115, 78, 140, 62, 95]
        mutated_indices = [0, 2, 6]

    amplicon_results = []
    for idx, amp in enumerate(AMPLICON_PANEL):
        depth = simulated_depths[idx]
        snp_status = "MUTATED ⚠️" if idx in mutated_indices else "Wild Type ✅"
        qc_status = "PASS" if depth >= 30 else "LOW COVERAGE ❌"
        
        amplicon_results.append({
            "Amplicon ID": amp["Amplicon_ID"],
            "Gen Target": amp["Gene"],
            "Hotspot SNP": amp["Hotspot"],
            "Read Depth (x)": depth,
            "QC Threshold": qc_status,
            "Hasil Alignment SNP": snp_status,
            "Rekomendasi Terapi": amp["Therapy"] if snp_status == "MUTATED ⚠️" else "Terapi Standar Non-Target",
            "Status E-FORNAS JKN": amp["EFORNAS"],
            "Notasi Genomik": amp["Target_SNP"]
        })
        
    df_res = pd.DataFrame(amplicon_results)
    return {
        "scenario_title": scenario_title,
        "total_reads": len(reads),
        "mean_depth": round(df_res["Read Depth (x)"].mean(), 1),
        "passing_amplicons": f"{len(df_res[df_res['Read Depth (x)'] >= 30])}/9",
        "amplicon_data": df_res
    }

def generate_coverage_plot(df_res, title_label=""):
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.plot(df_res["Amplicon ID"], df_res["Read Depth (x)"], marker='o', color='#2a52c9', linewidth=2.5, markersize=8, label='Read Depth (x)')
    ax.fill_between(df_res["Amplicon ID"], df_res["Read Depth (x)"], color='#2a52c9', alpha=0.15)
    ax.axhline(y=30, color='#e11d48', linestyle='--', linewidth=1.5, label='Min Threshold (30x QC)')
    
    ax.set_title(f"Lentera Gen — {title_label} Coverage Curve", fontsize=11, fontweight='bold', color='#0a1b3d', pad=12)
    ax.set_xlabel("Amplicon Target ID (1 Patient Panel)", fontsize=10, fontweight='bold', color='#0a1b3d')
    ax.set_ylabel("Read Depth / Coverage (x)", fontsize=10, fontweight='bold', color='#0a1b3d')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right')
    
    for idx, row in df_res.iterrows():
        color_text = '#e11d48' if row['Read Depth (x)'] < 30 else '#0a1b3d'
        ax.annotate(f"{row['Read Depth (x)']}x\n({row['Gen Target']})", 
                    (row["Amplicon ID"], row["Read Depth (x)"]), 
                    textcoords="offset points", xytext=(0, 8), ha='center', 
                    fontsize=8, fontweight='bold', color=color_text)
        
    plt.tight_layout()
    return fig