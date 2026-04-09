import streamlit as st
import primer3
import pandas as pd
from Bio import SeqIO  

# ==========================================
# 0. 网页全局配置 & 会话状态初始化
# ==========================================
st.set_page_config(
    page_title="引物设计平台 | SmartPrimer",
    page_icon="🧬",
    layout="wide", 
    initial_sidebar_state="expanded"
)

# 初始化引物购物车 (用于跨页面保存引物)
if 'primer_cart' not in st.session_state:
    st.session_state.primer_cart = []

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
    st.markdown("## ⚙️ 控制面板")
    
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
    st.markdown("#### 2. 参数配置")
    
    if is_qpcr:
        target_tm = st.slider("🎯 目标 Tm (°C)", 55.0, 65.0, 60.0, 0.5)
        amp_range = st.slider("📏 扩增子长度 (bp)", 50, 250, (70, 150), 10)
    else:
        if is_taichi: fragment_count = st.selectbox("🧩 总组装片段数 (含载体):", list(range(2, 17)), index=1)
        elif is_gibson: fragment_count = st.selectbox("🧩 总组装片段数 (含载体):", list(range(2, 7)), index=1)
        else: fragment_count = st.selectbox("🧩 总拼接片段数:", list(range(2, 7)), index=0)
        
        target_tm = st.slider("🎯 结合区目标 Tm (°C)", 50.0, 70.0, 65.0, 0.5)
        homology_len = st.number_input("🔗 同源臂长度 (bp)", 15, 60, 25)
            
    st.markdown("---")
    st.markdown("#### 3. 高级安全检查")
    do_enz_scan = st.checkbox("🧪 开启限制性内切酶扫描", value=True, help="扫描引物是否携带 EcoRI, BamHI 等 10 种常用酶切位点")

    st.markdown("<br><br><br>", unsafe_allow_html=True)
    st.caption("Powered by Primer3 Engine & AI")


# --- 主屏幕 (Main Content)：序列输入与结果展示 ---
st.title("🧬 智能引物设计平台")
st.markdown(f"**当前执行模式:** `< {method_choice} >`")
st.markdown("---")

st.markdown("### 📥 序列输入")
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


# --- 4. 执行计算与加入购物车 ---
st.markdown("<br>", unsafe_allow_html=True)
if st.button("🚀 开始设计 (并加入购物车)"):
    st.markdown("### 📊 本次引擎设计结果")
    if is_qpcr:
        if not gene_seq.strip(): st.error("⚠️ 请输入靶基因序列！")
        else:
            with st.spinner('🧬 正在扫描庞大的序列空间结构...'):
                results = design_qpcr_primers(gene_seq.replace(" ", "").replace("\n", "").upper(), target_tm, amp_range[0], amp_range[1], gene_name, do_enz_scan=do_enz_scan)
                if results:
                    st.success("🎉 设计完成！已为您筛选出表现最优的候选引物对，并加入购物车。")
                    
                    df = pd.DataFrame(results)
                    st.dataframe(df, use_container_width=True, hide_index=True)
                    
                    # ✅ 格式化并追加到全局购物车
                    for r in results:
                        st.session_state.primer_cart.append({"引物名称": r['正向引物'], "序列 (5'->3')": r["Fwd (5'->3')"], "备注": f"qPCR产物:{r['产物长']}bp"})
                        st.session_state.primer_cart.append({"引物名称": r['反向引物'], "序列 (5'->3')": r["Rev (5'->3')"], "备注": f"qPCR产物:{r['产物长']}bp"})
                else: st.error("❌ 未找到合适的引物对，请尝试放宽参数限制。")
    else:
        if any(not f['seq'] for f in all_fragments): st.error("⚠️ 请确保所有片段的序列框已填满！")
        else:
            with st.spinner(f'🔧 正在根据热力学规则规划 {current_method} 引物...'):
                results = design_assembly_primers(current_method, all_fragments, target_tm, homology_len, plasmid_name, do_enz_scan)
                df = pd.DataFrame(results)
                
                st.success("🎉 引物设计完成！已自动追加至下方购物车。")
                def highlight(val): return 'color: #00E5FF; font-weight: bold;' if isinstance(val, str) and '⚠️' in val else ''
                display_df = df.style.map(highlight, subset=['酶切警告']) if do_enz_scan else df
                st.dataframe(display_df, use_container_width=True, hide_index=True)
                
                # ✅ 格式化并追加到全局购物车
                for r in results:
                    st.session_state.primer_cart.append({"引物名称": r['引物名称'], "序列 (5'->3')": r["序列 (5'->3')"], "备注": r['备注']})


# ==========================================
# 5. 🛒 终极批量购物车面板
# ==========================================
st.markdown("<br><hr><br>", unsafe_allow_html=True)
st.markdown(f"### 🛒 引物批量暂存车 (当前共有 {len(st.session_state.primer_cart)} 条引物)")

if st.session_state.primer_cart:
    # 💡 新增：操作指南提示
    st.info("💡 **高阶操作指南**：\n"
            "1. **精准删除**：把鼠标移到表格最左侧的序号上，勾选对应行，点击右上角出现的『🗑️ 垃圾桶』即可删除个别引物。\n"
            "2. **直接修改**：双击任意单元格，可以像在 Excel 里一样，直接修改引物名称或备注信息！")
    
    cart_df = pd.DataFrame(st.session_state.primer_cart)
    
    # ✅ 核心升级：将不可交互的 dataframe 替换为强大的 data_editor
    edited_df = st.data_editor(
        cart_df,
        use_container_width=True,
        num_rows="dynamic", # 开启魔法参数：允许动态增加和删除行
        key="primer_editor"
    )
    
    # 实时将表格里的删除和修改动作，同步覆盖回全局购物车
    st.session_state.primer_cart = edited_df.to_dict('records')
    
    # 构建兼容 SnapGene 的批量 txt 格式 (使用编辑后的最新数据)
    txt_lines = [f"{row['引物名称']}\t{row['序列 (5\'->3\')']}\t{row['备注']}" for _, row in edited_df.iterrows()]
    txt_content = "\n".join(txt_lines)
    
    col_dl, col_clear, _ = st.columns([2, 1, 3])
    with col_dl:
        st.download_button(
            label="📥 批量导出全部引物 (SnapGene 格式)", 
            data=txt_content.encode('utf-8'), 
            file_name="Batch_Primers_Order.txt", 
            mime="text/plain"
        )
    with col_clear:
        if st.button("🚨 一键清空购物车"):
            st.session_state.primer_cart = []
            st.rerun()
else:
    st.info("💡 您的购物车空空如也，赶快去上面设计一些引物吧！不论切换什么质粒或模式，每次计算的结果都会在此自动累加。")
