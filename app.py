import streamlit as st
import primer3
import pandas as pd
from Bio import SeqIO  

# ==========================================
# 0. 网页全局配置 (必须放在第一行)
# ==========================================
st.set_page_config(
    page_title="智能引物设计平台 | SmartPrimer",
    page_icon="🧬",
    layout="wide", 
    initial_sidebar_state="expanded"
)

# --- 全新暗黑科技风 (Cyber-Genetics) CSS 注入 ---
st.markdown("""
<style>
    .stApp { background-color: #0E1117; }
    
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #161B22 0%, #0D1117 100%);
        border-right: 1px solid #30363D;
    }

    .stButton>button {
        width: 100%; border-radius: 8px; height: 50px;
        font-size: 18px; font-weight: 600; letter-spacing: 0.5px;
        background: linear-gradient(135deg, #6B21A8 0%, #00E5FF 100%);
        color: #FFFFFF; border: none;
        box-shadow: 0 4px 15px 0 rgba(0, 229, 255, 0.2);
        transition: all 0.3s ease 0s;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(0, 229, 255, 0.4);
        background: linear-gradient(135deg, #581C87 0%, #00B8D4 100%);
    }

    .stTextInput>div>div>input, .stSelectbox>div>div>div, .stTextArea>div>div>textarea {
        border-radius: 6px; border: 1px solid #30363D;
        background-color: #0D1117; color: #E2E8F0;
        transition: all 0.2s ease-in-out;
    }
    .stTextInput>div>div>input:focus, .stSelectbox>div>div>div:focus, .stTextArea>div>div>textarea:focus {
        border-color: #00E5FF; background-color: #161B22;
        box-shadow: 0 0 0 2px rgba(0, 229, 255, 0.25);
    }

    .stAlert {
        border-radius: 6px; border: none; border-left: 4px solid;
        background-color: #161B22; box-shadow: 0 1px 3px rgba(0,0,0,0.5);
    }

    /* 隐藏 Streamlit 默认菜单和底部水印 (保留 header 保证侧边栏按钮可用) */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)


# ==========================================
# 核心数据库：常用限制性内切酶
# ==========================================
COMMON_ENZYMES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
    "XhoI": "CTCGAG", "NdeI": "CATATG", "NotI": "GCGGCCGC",
    "SacI": "GAGCTC", "KpnI": "GGTACC", "SalI": "GTCGAC", "SpeI": "ACTAGT"
}

# ==========================================
# 核心算法一：文件解析引擎 (支持 SnapGene .dna)
# ==========================================
def parse_sequence_file(uploaded_file):
    sequences = []
    filename = uploaded_file.name.lower()
    
    if filename.endswith(".dna"):
        try:
            uploaded_file.seek(0)
            record = SeqIO.read(uploaded_file, "snapgene")
            seq_name = record.name if record.name and record.name != "<unknown name>" else filename.split('.')[0]
            sequences.append({"name": f"[完整] {seq_name}", "seq": str(record.seq).upper()})
            
            for feature in record.features:
                if feature.type in ["CDS", "promoter", "terminator", "misc_feature", "gene"]:
                    feature_name = feature.qualifiers.get('label', [feature.type])[0]
                    feature_seq = str(feature.extract(record.seq)).upper()
                    sequences.append({"name": f"[元件] {feature_name}", "seq": feature_seq})
        except Exception:
            pass 
    else:
        uploaded_file.seek(0)
        file_content = uploaded_file.getvalue().decode("utf-8")
        lines = file_content.splitlines()
        
        if not lines: return sequences
        if lines[0].strip().startswith(">"):
            curr_name, curr_seq = "", []
            for line in lines:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if curr_name: sequences.append({"name": curr_name, "seq": "".join(curr_seq).upper()})
                    curr_name, curr_seq = line[1:].strip(), []
                else:
                    curr_seq.append(line)
            if curr_name: sequences.append({"name": curr_name, "seq": "".join(curr_seq).upper()})
        else:
            raw_seq = "".join([line.strip() for line in lines]).upper()
            if raw_seq: sequences.append({"name": filename.split('.')[0], "seq": raw_seq})
            
    return sequences

# ==========================================
# 核心算法二：引物设计辅助引擎
# ==========================================
def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def check_restriction_sites(sequence):
    found_enzymes = [e for e, site in COMMON_ENZYMES.items() if site in sequence]
    return "⚠️ " + ", ".join(found_enzymes) if found_enzymes else "✅ 安全"

def get_binding_sequence(seq, target_tm, from_5_prime=True):
    for length in range(15, 46):
        binding_seq = seq[:length] if from_5_prime else seq[-length:]
        if primer3.calc_tm(binding_seq) >= target_tm: return binding_seq, primer3.calc_tm(binding_seq)
    fallback_seq = seq[:45] if from_5_prime else seq[-45:]
    return fallback_seq, primer3.calc_tm(fallback_seq)

def design_assembly_primers(method, all_fragments, target_tm, homology_len, plasmid_name="", do_enz_scan=True):
    primers_list = []
    N = len(all_fragments)
    is_circular = method in ["Gibson", "Tai Chi"] 

    for i in range(N):
        curr_name = all_fragments[i]['name']
        curr_seq = all_fragments[i]['seq']
        curr_template = all_fragments[i].get('template', curr_name)
        kb_len = f"{len(curr_seq)/1000:.1f}k"
        
        if is_circular:
            fwd_primer_name = f"{plasmid_name}-{N}-{i+1}-{kb_len}-{curr_template}-F"
            rev_primer_name = f"{plasmid_name}-{N}-{i+1}-{kb_len}-{curr_template}-R"
        else:
            fwd_primer_name = f"{curr_name}-F"
            rev_primer_name = f"{curr_name}-R"
        
        # Fwd Primer
        bind_fwd, tm_fwd = get_binding_sequence(curr_seq, target_tm, from_5_prime=True)
        if is_circular:
            prev_idx = (i - 1) % N
            homology_seq = all_fragments[prev_idx]['seq'][-homology_len:]
            primer_fwd = homology_seq + bind_fwd
            note_fwd = f"连 [{all_fragments[prev_idx]['name']}] 3'端"
        else:
            if i == 0:
                primer_fwd, note_fwd = bind_fwd, "Overlap: 最外侧 Fwd"
            else:
                homology_seq = all_fragments[i-1]['seq'][-homology_len:]
                primer_fwd = homology_seq + bind_fwd
                note_fwd = f"连 [{all_fragments[i-1]['name']}] 3'端"
            
        fwd_data = {
            "序号": len(primers_list) + 1, "引物名称": fwd_primer_name, "序列 (5'->3')": primer_fwd, 
            "长度": len(primer_fwd), "Tm": round(tm_fwd, 2), "备注": note_fwd
        }
        if do_enz_scan: fwd_data["酶切警告"] = check_restriction_sites(primer_fwd)
        primers_list.append(fwd_data)

        # Rev Primer
        bind_rev, tm_rev = get_binding_sequence(curr_seq, target_tm, from_5_prime=False)
        if is_circular:
            next_idx = (i + 1) % N
            homology_seq = get_reverse_complement(all_fragments[next_idx]['seq'][:homology_len])
            primer_rev = homology_seq + get_reverse_complement(bind_rev)
            note_rev = f"连 [{all_fragments[next_idx]['name']}] 5'端"
        else:
            if i == N - 1:
                primer_rev, note_rev = get_reverse_complement(bind_rev), "Overlap: 最外侧 Rev"
            else:
                homology_seq = get_reverse_complement(all_fragments[i+1]['seq'][:homology_len])
                primer_rev = homology_seq + get_reverse_complement(bind_rev)
                note_rev = f"连 [{all_fragments[i+1]['name']}] 5'端"
            
        rev_data = {
            "序号": len(primers_list) + 1, "引物名称": rev_primer_name, "序列 (5'->3')": primer_rev, 
            "长度": len(primer_rev), "Tm": round(tm_rev, 2), "备注": note_rev
        }
        if do_enz_scan: rev_data["酶切警告"] = check_restriction_sites(primer_rev)
        primers_list.append(rev_data)
        
    return primers_list

def design_qpcr_primers(target_seq, target_tm, min_amp, max_amp, gene_name, max_results=3, do_enz_scan=True):
    results = []
    seq_len = len(target_seq)
    for i in range(seq_len - max_amp):
        for f_len in range(18, 23):
            f_seq = target_seq[i : i+f_len]
            f_tm = primer3.calc_tm(f_seq)
            if abs(f_tm - target_tm) > 1.5 or primer3.calc_hairpin(f_seq).structure_found or primer3.calc_homodimer(f_seq).structure_found: continue
            
            for amp_len in range(min_amp, max_amp + 1):
                r_start = i + amp_len
                if r_start > seq_len: break
                for r_len in range(18, 23):
                    r_seq = get_reverse_complement(target_seq[r_start-r_len : r_start])
                    r_tm = primer3.calc_tm(r_seq)
                    if abs(r_tm - target_tm) > 1.5 or abs(f_tm - r_tm) > 1.0: continue
                    if primer3.calc_hairpin(r_seq).structure_found or primer3.calc_homodimer(r_seq).structure_found or primer3.calc_heterodimer(f_seq, r_seq).structure_found: continue
                    
                    pair_idx = len(results) + 1
                    pair_data = {
                        "序号": pair_idx, "方案": f"Pair {pair_idx}", 
                        "正向引物": f"{gene_name}-{pair_idx}-F", "Fwd (5'->3')": f_seq, "Fwd Tm": round(f_tm, 2),
                        "反向引物": f"{gene_name}-{pair_idx}-R", "Rev (5'->3')": r_seq, "Rev Tm": round(r_tm, 2), 
                        "产物长": amp_len, "Tm 差": round(abs(f_tm - r_tm), 2)
                    }
                    if do_enz_scan:
                        pair_data["Fwd 酶切"] = check_restriction_sites(f_seq)
                        pair_data["Rev 酶切"] = check_restriction_sites(r_seq)
                    results.append(pair_data)
                    if len(results) >= max_results: return results
    return results


# ==========================================
# 网页界面 (UI 设计区)
# ==========================================

# --- 侧边栏 (Sidebar)：所有设置与控制 ---
with st.sidebar:
    st.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=60)
    st.markdown("## ⚙️ 核心控制面板")
    
    st.markdown("#### 1. 实验方案选择")
    method_choice = st.radio(
        "请选择:", 
        ["Tai Chi Assembly (太极组装)", "Gibson Assembly (多片段闭环)", "Overlap PCR (线性拼接)", "qPCR (定量 PCR)"]
    )
    
    is_taichi = "Tai Chi" in method_choice
    is_gibson = "Gibson" in method_choice
    is_overlap = "Overlap" in method_choice
    is_qpcr = "qPCR" in method_choice
    needs_vector = is_taichi or is_gibson 

    if is_taichi: current_method = "Tai Chi"
    elif is_gibson: current_method = "Gibson"
    elif is_overlap: current_method = "Overlap"
    else: current_method = "qPCR"

    st.markdown("---")
    st.markdown("#### 2. 热力学参数配置")
    
    if is_qpcr:
        target_tm = st.slider("🎯 目标 Tm (°C)", 55.0, 65.0, 60.0, 0.5)
        amp_range = st.slider("📏 扩增子长度 (bp)", 50, 250, (70, 150), 10)
    else:
        if is_taichi: fragment_count = st.selectbox("🧩 总组装片段数 (含载体):", list(range(2, 17)), index=1)
        elif is_gibson: fragment_count = st.selectbox("🧩 总组装片段数 (含载体):", list(range(2, 7)), index=1)
        else: fragment_count = st.selectbox("🧩 总拼接片段数:", list(range(2, 7)), index=0)
        
        target_tm = st.slider("🎯 结合区目标 Tm (°C)", 50.0, 70.0, 60.0, 0.5)
        homology_len = st.number_input("🔗 同源臂长度 (bp)", 15, 60, 25)
            
    st.markdown("---")
    st.markdown("#### 3. 高级安全检查")
    do_enz_scan = st.checkbox("🧪 开启限制性内切酶扫描", value=True, help="扫描引物是否携带 EcoRI, BamHI 等 10 种常用酶切位点")

    st.markdown("<br><br><br>", unsafe_allow_html=True)
    st.caption("Powered by Primer3 Engine & AI")


# --- 主屏幕 (Main Content)：序列输入与结果展示 ---
st.title("🧬 智能核酸引物设计平台")
st.markdown(f"**当前执行模式:** `< {method_choice} >`")
st.markdown("---")

st.markdown("### 📥 序列输入舱")
uploaded_file = st.file_uploader("📂 支持拖拽 SnapGene .dna, .fasta 或 .txt，实现自动填表与元件提取", type=["fasta", "fas", "txt", "seq", "dna"])

imported_seqs = []
if uploaded_file is not None:
    imported_seqs = parse_sequence_file(uploaded_file)
    if imported_seqs: st.toast(f"✅ 成功提取 {len(imported_seqs)} 条序列！")

plasmid_name, gene_name, all_fragments = "", "", []

if is_qpcr:
    gene_name = st.text_input("🏷️ 靶基因名称", value=imported_seqs[0]["name"].replace("[完整] ", "").replace("[元件] ", "") if imported_seqs else "Target_Gene")
    gene_seq = st.text_area("🧬 靶基因完整序列 (5' -> 3')", value=imported_seqs[0]["seq"] if imported_seqs else "", height=200)

else:
    if needs_vector:
        plasmid_name = st.text_input(f"🏆 最终构建的质粒名称 (用于订单命名):", value="pNew_Plasmid")
        
        st.info("💡 **片段 1 (载体骨架)** 将被视为组装的基石，平台会自动为其设计线性化扩增引物。")
        v_col1, v_col2 = st.columns([1, 3])
        with v_col1: 
            v_name = st.text_input("载体命名", value=imported_seqs[0]["name"].replace("[完整] ", "").replace("[元件] ", "") if imported_seqs else "Vector")
            v_temp = st.text_input("🧪 扩增用 PCR 模板", value=v_name, key="v_temp")
        with v_col2: 
            v_seq = st.text_area("载体序列 (5' -> 3')", value=imported_seqs[0]["seq"] if imported_seqs else "", height=130)
        
        all_fragments.append({"name": v_name.strip(), "seq": v_seq.replace(" ", "").replace("\n", "").upper(), "template": v_temp.strip()})

        st.markdown("#### 🧩 插入片段序列")
        for i in range(1, fragment_count):
            d_name = imported_seqs[i]["name"].replace("[完整] ", "").replace("[元件] ", "") if i < len(imported_seqs) else f"Insert_{i}"
            d_seq = imported_seqs[i]["seq"] if i < len(imported_seqs) else ""
            f_col1, f_col2 = st.columns([1, 3])
            with f_col1: 
                f_name = st.text_input(f"片段 {i+1} 命名", value=d_name, key=f"fn_{i}")
                f_temp = st.text_input("🧪 扩增用 PCR 模板", value=d_name, key=f"ft_{i}")
            with f_col2: 
                f_seq = st.text_area(f"片段 {i+1} 序列 (5' -> 3')", value=d_seq, height=130, key=f"fs_{i}")
            
            all_fragments.append({"name": f_name.strip(), "seq": f_seq.replace(" ", "").replace("\n", "").upper(), "template": f_temp.strip()})
    else:
        st.markdown("#### 🧩 线性拼接片段序列")
        for i in range(fragment_count):
            d_name = imported_seqs[i]["name"].replace("[完整] ", "").replace("[元件] ", "") if i < len(imported_seqs) else f"Fragment_{i+1}"
            d_seq = imported_seqs[i]["seq"] if i < len(imported_seqs) else ""
            f_col1, f_col2 = st.columns([1, 3])
            with f_col1: 
                f_name = st.text_input(f"片段 {i+1} 命名", value=d_name, key=f"fn_{i}")
                f_temp = st.text_input("🧪 扩增用 PCR 模板", value=d_name, key=f"ft_{i}")
            with f_col2: 
                f_seq = st.text_area(f"片段 {i+1} 序列 (5' -> 3')", value=d_seq, height=130, key=f"fs_{i}")
                
            all_fragments.append({"name": f_name.strip(), "seq": f_seq.replace(" ", "").replace("\n", "").upper(), "template": f_temp.strip()})


# --- 4. 执行计算 ---
st.markdown("<br>", unsafe_allow_html=True)
if st.button("🚀 启动 AI 引擎进行设计"):
    st.markdown("### 📊 引擎设计结果")
    if is_qpcr:
        if not gene_seq.strip(): st.error("⚠️ 请输入靶基因序列！")
        else:
            with st.spinner('🧬 正在扫描庞大的序列空间结构...'):
                results = design_qpcr_primers(gene_seq.replace(" ", "").replace("\n", "").upper(), target_tm, amp_range[0], amp_range[1], gene_name, do_enz_scan=do_enz_scan)
                if results:
                    st.success("🎉 设计完成！已为您筛选出表现最优的候选引物对。")
                    
                    df = pd.DataFrame(results)
                    # ✅ 这里加上 hide_index=True 隐藏 Pandas 自带序号
                    st.dataframe(df, use_container_width=True, hide_index=True)
                    
                    txt_lines = []
                    for r in results:
                        txt_lines.append(f"{r['正向引物']}\t{r['Fwd (5\'->3\')']}\tqPCR产物:{r['产物长']}bp")
                        txt_lines.append(f"{r['反向引物']}\t{r['Rev (5\'->3\')']}\tqPCR产物:{r['产物长']}bp")
                    txt_content = "\n".join(txt_lines)
                    
                    st.download_button(
                        label="📥 下载 SnapGene 引物导入文件 (.txt)", 
                        data=txt_content.encode('utf-8'), 
                        file_name=f"{gene_name}_qPCR_Primers.txt", 
                        mime="text/plain"
                    )
                else: st.error("❌ 未找到合适的引物对，请尝试放宽参数限制。")
    else:
        if any(not f['seq'] for f in all_fragments): st.error("⚠️ 请确保所有片段的序列框已填满！")
        else:
            with st.spinner(f'🔧 正在根据热力学规则规划 {current_method} 引物...'):
                results = design_assembly_primers(current_method, all_fragments, target_tm, homology_len, plasmid_name, do_enz_scan)
                df = pd.DataFrame(results)
                
                st.success("🎉 所有引物设计完成！请检查下方的酶切位点警告状态。")
                def highlight(val): return 'color: #00E5FF; font-weight: bold;' if isinstance(val, str) and '⚠️' in val else ''
                
                display_df = df.style.map(highlight, subset=['酶切警告']) if do_enz_scan else df
                # ✅ 这里加上 hide_index=True 隐藏 Pandas 自带序号
                st.dataframe(display_df, use_container_width=True, hide_index=True)
                
                txt_lines = []
                for r in results:
                    txt_lines.append(f"{r['引物名称']}\t{r['序列 (5\'->3\')']}\t{r['备注']}")
                txt_content = "\n".join(txt_lines)
                
                fname = f"{plasmid_name}_{current_method.replace(' ', '')}.txt" if needs_vector else f"Overlap_{fragment_count}Frags.txt"
                st.download_button(
                    label="📥 下载 SnapGene 引物导入文件 (.txt)", 
                    data=txt_content.encode('utf-8'), 
                    file_name=fname, 
                    mime="text/plain"
                )
