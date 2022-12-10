import streamlit as st
import pandas as pd
import regex as re
import time
import io
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import gc_fraction


def app():
    st.markdown("# amiR Design Tool")
    st.markdown("You can use this app to design artificial miRNA to knock down your favourite gene.")
    st.markdown("**Paste your target CDS in the box below, select a Database to BLAST against, then press START**.")
    st.markdown("<style>div.stButton > button:first-child {background-color: #0099ff;color:#ffffff;}</style>", unsafe_allow_html=True)    
    
    def validate(seq):
        leftover = set(seq) - set("ATGCN")
        return not leftover

    def search_filter_oligos(seq):
        oligos = re.findall(r"[C|G].[^A]......[A|T]...........[A|T]", seq, overlapped=True, flags=re.DOTALL) 
        oligos = list(filter(lambda seq: 'AAAAA' not in seq and 'CCCC' not in seq and 'GGGG' not in seq and 'TTTTT' not in seq, oligos)) 
        return [oligo for oligo in oligos if gc_fraction(oligo[0:9]) > 0.5 and gc_fraction(oligo[9:22]) < 0.5] 

    def clear_form():
        st.session_state["user_input"] = " "
        st.session_state["user_select"] = "Mouse"
    

    user_input_seq = st.text_area(label=" ", key="user_input")
    user_input_seq = "".join(letter.rstrip() for letter in user_input_seq).upper()
    
    select_db_options = {"mus_mod_database": "Mouse","h_mod_database": "Human"}
    select_db = st.radio("Please select the preferred Database",options=select_db_options.values(), key="user_select")
    select_db = next((key for key, value in select_db_options.items() if value == select_db), None)

    col1, col2 = st.columns([0.42,4])
    with col1:
        start = st.button("Start")
    with col2:
        st.button("Reset Input",on_click=clear_form)
        
        
    if start:
        if "".__eq__(user_input_seq) or not validate(user_input_seq):
            st.markdown("<font color=red>**Please insert a valid CDS sequence**</font>", unsafe_allow_html=True)
        else:
            start= time.time()

            filtered_oligos = search_filter_oligos(user_input_seq)
            oligos = {index:[filtered_oligos[index],filtered_oligos[index][14:21]] for index in range(0,len(filtered_oligos))} ## SEED PONTOSAN HOL VAN? [14:21]??
            oligodf = pd.DataFrame.from_dict(oligos, orient="index", columns=["oligo sequence","seed sequence"])
            oligodf['id'] = oligodf.index

            fasta = "\n".join([">" + str(row[0]) + "\n" + row[1] + "\n" for row in oligodf.itertuples()])

            blast_results = NcbiblastnCommandline(db = select_db, word_size=7, dust ="no", evalue=10, num_threads=2, outfmt ="6 qseqid sseqid bitscore sseq") 
            out, err = blast_results(stdin=fasta)
            
            blast_results_df = pd.read_csv(io.StringIO(out), sep="\t", header=None, names=["id", "sseqid", "bitscore", "sseq"])
            
            blast_results_df = pd.merge(blast_results_df.astype(str), oligodf.astype(str), on="id")
            blast_results_df.replace(to_replace=r"'", value=" ", regex=True, inplace=True)
            blast_results_df.replace(to_replace=r"\.\d", value=" ", regex=True, inplace=True)
            blast_results_df[['sseqid','gene description']] = blast_results_df["sseqid"].str.split(" ", n = 1, expand=True)

            blast_results_df = blast_results_df[:].loc[blast_results_df.apply(lambda x: x["seed sequence"] in x["sseq"], axis=1)]

            tophits = pd.DataFrame(blast_results_df.groupby("id").size().nsmallest(5, keep="all"))
            tophitsout = pd.DataFrame()

            for row in tophits.iterrows():

                tophitsout = pd.concat([tophitsout,(blast_results_df[:].loc[lambda df: df["id"] == row[0]])],ignore_index=True)

            st.write((time.time() - start))
            st.write(tophitsout)

                    
if __name__ == "__main__":
    app()