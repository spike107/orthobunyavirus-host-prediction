import os
import sys
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import re
import openpyxl
from difflib import get_close_matches
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  
import matplotlib.patches as patches  
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
import json
from sklearn.metrics import classification_report, accuracy_score, f1_score, recall_score, roc_auc_score
import umap

try:
    import iFeatureOmegaCLI 
    HAS_IFEATURE = True
except Exception:
    HAS_IFEATURE = False

import warnings
warnings.filterwarnings('ignore')
from datetime import datetime
import argparse
from pathlib import Path
import shutil

# ========== Tee class: write to file and console simultaneously ==========
class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()

# =====================================================
# Phase 1: Build the main label table (metadata + final label)
# =====================================================

def normalize_name(s):
    """Normalise virus names: 
    lowercase, remove 'virus'/'orthobunyavirus' suffixes, 
    strip non-alphanumerics."""
    s = str(s).lower()
    s = s.replace('orthobunyavirus', '').replace('virus', '')
    s = re.sub(r'[^a-z0-9 ]+', '', s)
    return s.strip()

def load_references():
    """Load external reference resources."""
    # ICTV species list.
    ictv_species = pd.read_csv('ICTV_2024_MSL.csv', encoding='utf-8')
    ictv_species_set = set(ictv_species['Species'].astype(str).str.strip())
    
    # ICTV rename table.
    taxa_renamed = pd.read_csv('ICTV_2024_Taxa Renamed or Abolished.csv', encoding='utf-8')
    taxa_renamed = taxa_renamed.dropna(subset=['Old Name', 'New Name'])
    old2new = dict(zip(taxa_renamed['Old Name'].str.strip(), taxa_renamed['New Name'].str.strip()))
    
    # Virus-Host DB list of human viruses.
    with open('human_virus_DB_species_clean.txt', encoding='utf-8') as f:
        virus_hostdb_human_set = set(normalize_name(line.strip()) for line in f if line.strip())
    
    return ictv_species_set, old2new, virus_hostdb_human_set

def load_ictv_orthobunyavirus_mapping(csv_file):
    """Load ICTV Orthobunyavirus mapping: 
    new Latin name → legacy virus name."""
    ictv = pd.read_csv(csv_file)
    latin2virus = {
        normalize_name(row['Latin_name']): row['Virus_name'].strip()
        for _, row in ictv.iterrows()
        if pd.notna(row['Latin_name']) and pd.notna(row['Virus_name'])
    }
    return latin2virus

def load_orthobunyavirus_taxonomy(file_path):
    """Parse the Orthobunyavirus_taxonomy.xlsx structure."""
    wb = openpyxl.load_workbook(file_path)
    ws = wb.active
    species_map = {}
    group = None
    
    for row in ws.iter_rows(min_row=2):
        cell = row[0].value
        if not cell:
            continue
        if row[0].font.italic:
            group = cell.strip()
            continue
        strain = cell.strip()
        strain_main = strain.split('-')[0].strip()
        strain_norm = normalize_name(strain_main)
        species_map[strain_norm] = {'group': group if group else ''}
    
    return species_map

def standardize_species_with_manual(s, manual_species_dict, old2new, ictv_species_set, latin2virus):
    """Standardise species name; 
    return (standardised name, evidence source)."""
    s_main = str(s).split('-')[0].strip()
    s_norm = normalize_name(s_main)
    # Priority matching order.
    # 1) ICTV new (Latin) name → legacy virus name.
    if s_norm in latin2virus:
        return latin2virus[s_norm], 'ictv_new2old'
    # 2) Manual mapping table: exact match.
    if s_norm in manual_species_dict:
        return manual_species_dict[s_norm]['group'], 'manual_exact'
    # 3) Containment (substring) match.
    for k in manual_species_dict.keys():
        if len(s_norm) >= 6 and (s_norm in k or k in s_norm):
            if len(s_norm) < len(k)*1.8 and len(k) < len(s_norm)*1.8:
                return manual_species_dict[k]['group'], 'manual_contain'
    # 4) Fuzzy match.
    close = get_close_matches(s_norm, manual_species_dict.keys(), n=1, cutoff=0.8)
    if close:
        return manual_species_dict[close[0]]['group'], 'manual_fuzzy'
    # 5) ICTV old → new name.
    old2new_norm = {normalize_name(k): v for k, v in old2new.items()}
    if s_norm in old2new_norm:
        return old2new_norm[s_norm], 'ictv_old2new'
    # 6) In ICTV species set.
    ictv_species_set_norm = set(normalize_name(x) for x in ictv_species_set)
    if s_norm in ictv_species_set_norm:
        matched = [x for x in ictv_species_set if normalize_name(x) == s_norm]
        return matched[0] if matched else s, 'ictv'
    # 7) Unmatched.
    return s, 'unmatched'

def assign_label(species_std, virus_hostdb_human_set):
    """Assign human/nonhuman label by standardised species name."""
    s_norm = normalize_name(species_std)
    return 'human' if s_norm in virus_hostdb_human_set else 'nonhuman'

def parse_nucleotide_header(record):
    """Extract information from nucleotide FASTA headers."""
    desc = record.description
    # Extract segment.
    seg = desc.strip().split('|')[-1].strip()
    segment = seg if seg in ['S', 'M', 'L'] else 'Unknown'
    # Extract species.
    def extract_main_species(s):
        words = s.strip().split()
        return ' '.join(words[:2]) if len(words) >= 2 else s.strip()
    try:
        candidates = desc.split('|')
        if len(candidates) > 3:
            field = candidates[3].strip()
            if "virus" in field.lower() or re.match(r"^[A-Z][a-z]+ [a-z]+", field):
                species = extract_main_species(field)
            else:
                species = "Unknown"
        else:
            for field in candidates:
                if "virus" in field.lower() or re.match(r"^[A-Z][a-z]+ [a-z]+", field):
                    species = extract_main_species(field)
                    break
            else:
                species = "Unknown"
    except:
        species = "Unknown"
    # Extract completeness.
    desc_lower = desc.lower()
    completeness_raw = 'other'
    for kw in ['complete genome', 'complete sequence', 'complete cds', 'partial cds', 'partial sequence', 'gene']:
        if kw in desc_lower:
            completeness_raw = kw
            break
    # Normalise completeness labels.
    completeness_map = {
        'complete sequence': 'complete_sequence',
        'complete genome': 'complete_sequence',
        'complete cds': 'complete_cds',
        'partial cds': 'partial_cds',
        'partial sequence': 'partial_sequence',
        'gene': 'gene',
        'other': 'other'
    }
    completeness = completeness_map.get(completeness_raw, 'other')
    return {
        'id': record.id,
        'species': species,
        'segment': segment,
        'completeness': completeness
    }

def build_meta_main(human_fasta, nonhuman_fasta, output_csv='meta_main_post_relabel.csv'):
    print("\n=== Phase 1: Building Main Label Table (No Completeness Filter) ===")  
    # Load reference resources.
    ictv_species_set, old2new, virus_hostdb_human_set = load_references()
    manual_species_dict = load_orthobunyavirus_taxonomy("Orthobunyavirus_taxonomy.xlsx")
    latin2virus = load_ictv_orthobunyavirus_mapping('ICTV Orthobunyavirus.csv')
    all_records = []
    
    # Process human FASTA.
    print("Processing human FASTA...")
    for record in SeqIO.parse(human_fasta, 'fasta'):
        info = parse_nucleotide_header(record)
        info['source_file'] = 'human'
        info['species_std'] = info['species']  
        info['label_evidence'] = 'original_human'
        info['final_label'] = 'human'
        all_records.append(info)
    
    # Process nonhuman FASTA.
    print("Processing nonhuman FASTA...")
    for record in SeqIO.parse(nonhuman_fasta, 'fasta'):
        info = parse_nucleotide_header(record)
        info['source_file'] = 'nonhuman'
        
        # Standardise species names.
        species_std, evidence = standardize_species_with_manual(
            info['species'], manual_species_dict, old2new, ictv_species_set, latin2virus
        )
        info['species_std'] = species_std
        info['label_evidence'] = evidence
        
        # Assign labels.
        info['final_label'] = assign_label(species_std, virus_hostdb_human_set)
        all_records.append(info)
    
    # Create DataFrame.
    meta_main = pd.DataFrame(all_records)
    
    # Clean suspicious species records.
    print("Cleaning suspicious species...")
    suspicious_names = ['Botrytis cinerea', 'Centruroides sculpturatus', 'Diaphorina citri']
    removed = []
    for name in suspicious_names:
        mask = meta_main['species'].str.contains(name, case=False, na=False)
        if mask.any():
            print(f"  Removing {mask.sum()} records containing '{name}'")
            removed_records = meta_main[mask].copy()
            removed_records['removed_reason'] = name
            removed.append(removed_records)
            meta_main = meta_main[~mask]
    if removed:
        pd.concat(removed).to_csv("removed_suspicious_species.csv", index=False)
        print("✓ Removed suspicious species records saved: removed_suspicious_species.csv")
        
    # **Remove completeness filtering** - keep all completeness types.
    print("Keeping all sequences regardless of completeness...")
    print(f"Completeness distribution:")
    completeness_dist = meta_main['completeness'].value_counts()
    for comp_type, count in completeness_dist.items():
        print(f"  {comp_type}: {count}")

    # Fix apostrophe encoding issues.
    print("Fixing encoding issues...")
    for col in ['species', 'species_std']:
        meta_main[col] = meta_main[col].str.replace(''', "'").str.replace(''', "'")
        meta_main[col] = (meta_main[col]
                         .str.replace('鈥橮', "'", regex=False)  # M鈥橮oko → M'oko
                         .str.replace('谩', 'á', regex=False)   # Guajar谩 → Guajará  
                         .str.replace('艌', 'Ń', regex=False))  # 扭ahy艌a → ŃahyŃa    
    # Save main table. 
    output_csv = 'meta_main_post_relabel.csv'
    meta_main.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"✓ Main metadata table saved: {output_csv}")
    print(f"  Total records: {len(meta_main)}")
    print(f"  Human sequences: {(meta_main['final_label']=='human').sum()}")
    print(f"  Nonhuman sequences: {(meta_main['final_label']=='nonhuman').sum()}")

    # 1) Create the core species summary table (Species Summary).
    print("Generating master species summary table...")
    
    # Group by species.
    species_grouped = meta_main.groupby('species_std')
    
    # Build summary records.
    summary_data = []
    for species, group in species_grouped:
        human_count = (group['final_label'] == 'human').sum()
        nonhuman_count = (group['final_label'] == 'nonhuman').sum()
        
        # Decide species-level final label.
        species_final_label = 'human' if human_count > 0 else 'nonhuman'
        
        # Summarise label evidence and source files.
        # Consider only evidence that determines the final label.
        if species_final_label == 'human':
            relevant_evidence = group[group['final_label'] == 'human']['label_evidence'].unique()
        else:
            relevant_evidence = group['label_evidence'].unique() # For nonhuman, include all evidence.
            
        evidence_summary = ', '.join(sorted(list(relevant_evidence)))
        source_file_summary = ', '.join(sorted(list(group['source_file'].unique())))

        summary_data.append({
            'species_std': species,
            'species_final_label': species_final_label,
            'seq_count_human': human_count,
            'seq_count_nonhuman': nonhuman_count,
            'total_seqs': len(group),
            'label_evidence_summary': evidence_summary,
            'source_file_summary': source_file_summary
        })

    species_summary_df = pd.DataFrame(summary_data)
    species_summary_df = species_summary_df.sort_values(by=['species_final_label', 'total_seqs'], ascending=[True, False])
    
    # Save species summary table.
    species_summary_output = "species_summary_final.csv"
    species_summary_df.to_csv(species_summary_output, index=False, encoding='utf-8')
    print(f"✓ Master species summary saved: {species_summary_output}")
    print(f"  Total unique species: {len(species_summary_df)}")
    print(f"  Species classified as 'human': {(species_summary_df['species_final_label'] == 'human').sum()}")
    print(f"  Species classified as 'nonhuman': {(species_summary_df['species_final_label'] == 'nonhuman').sum()}")

    # 2) Generate relabelling statistics.
    print("Generating relabeling summary...")
    relabel_stats = meta_main.groupby(['source_file', 'final_label']).size().unstack(fill_value=0)
    relabel_stats.to_csv("relabel_summary.csv")
    print("✓ Relabeling summary saved: relabel_summary.csv")
    if 'human' in relabel_stats.columns and 'nonhuman' in relabel_stats.index:
        relabelled_count = relabel_stats.loc['nonhuman', 'human']
        print(f"  > {relabelled_count} sequences from nonhuman file were relabeled to 'human'.")

    # 3) Unmatched-species summar
    print("Processing unmatched species reports...")
    unmatched = meta_main[meta_main['label_evidence'] == 'unmatched']
    if not unmatched.empty:
        unmatched[['id', 'species']].to_csv("unmatched_species.csv", index=False)
        print(f"✓ Unmatched species list saved: unmatched_species.csv ({len(unmatched)} records)")
        # Count unmatched records by species.
        unmatched_species_counts = (
            unmatched.groupby('species_std')
            .size()
            .reset_index(name='count')
            .sort_values(by='count', ascending=False)
        )
        unmatched_species_counts.to_csv("unmatched_species_counts.csv", index=False)
        print(f"✓ Unmatched species count summary saved: unmatched_species_counts.csv ({len(unmatched_species_counts)} unique species)")
    else:
        print("✓ No unmatched species found.")

    return meta_main

# =====================================================
# Phase 2: Preprocess and align the three FASTA sets.
# =====================================================
def parse_cds_header(record):
    """Extract info from CDS FASTA headers — custom pipe-delimited format."""
    info = {
        'cds_id': record.id,
        'protein_id': None,
        'nucleotide_id': None,
        'nucleotide_id_clean': None,
        'completeness': None
    }
    # Extract nucleotide ID from record.id (before the colon).
    id_parts = record.id.split(':')
    if len(id_parts) >= 1:
        info['nucleotide_id'] = id_parts[0]
        info['nucleotide_id_clean'] = re.sub(r'\.\d+$', '', id_parts[0])
    # Parse the pipe-delimited description.
    # Format：>NC_077898.1:32..6787 |description|protein_id|genome_id|source|virus_name|species|segment|length|completeness|host
    desc_parts = record.description.split('|')
    if len(desc_parts) >= 3:
        # protein_id is field 3 (index 2).
        potential_prot_id = desc_parts[2].strip()
        if re.match(r'[A-Z]{2,3}_[0-9]{6,}\.\d+', potential_prot_id):
            info['protein_id'] = potential_prot_id
    if len(desc_parts) >= 10:
        # completeness is field 10 (index 9).
        completeness_field = desc_parts[9].strip().lower()
        if completeness_field == 'complete':
            info['completeness'] = 'complete_cds'
        elif completeness_field == 'partial':
            info['completeness'] = 'partial_cds'
        else:
            info['completeness'] = 'unknown'
    return info

def parse_protein_header(record):
    """Parse protein FASTA headers."""
    desc = record.description
    parts = desc.split('|')
    
    info = {
        'protein_id': record.id,
        'nucleotide_id': parts[4].strip() if len(parts) > 4 else None,
    }
    # Drop version suffix.
    if info['nucleotide_id']:
        info['nucleotide_id_clean'] = re.sub(r'\.\d+$', '', info['nucleotide_id'])
    
    return info

def label_cds_fasta(cds_fasta, meta_main, output_csv='cds_meta_labeled.csv'):
    # Keep only records whose labels are found in meta_main.
    # Use for ancillary description/metadata and reporting; no longer used as the base input for RSCU features.
    print("\n=== Processing CDS FASTA ===")
    # Prepare mapping tables.
    meta_main['id_clean'] = meta_main['id'].str.replace(r'\.\d+$', '', regex=True)
    id_to_label = dict(zip(meta_main['id_clean'], meta_main['final_label']))
    id_to_species_std = dict(zip(meta_main['id_clean'], meta_main['species_std']))
    id_to_segment = dict(zip(meta_main['id_clean'], meta_main['segment']))
    
    cds_records = []
    for record in SeqIO.parse(cds_fasta, 'fasta'):
        info = parse_cds_header(record)
        
        # Map labels via nucleotide_id.
        if info.get('nucleotide_id_clean'):
            info['species_std'] = id_to_species_std.get(info['nucleotide_id_clean'], 'Unknown')
            info['final_label'] = id_to_label.get(info['nucleotide_id_clean'], 'Unknown')
            info['segment'] = id_to_segment.get(info['nucleotide_id_clean'], 'Unknown')
        else:
            info['species_std'] = 'Unknown'
            info['final_label'] = 'Unknown'
            info['segment'] = 'Unknown'
        
        info['seq_length'] = len(record.seq)
        cds_records.append(info)
    
    cds_meta = pd.DataFrame(cds_records)
    cds_meta = cds_meta[cds_meta['final_label'] != 'Unknown']
    
    cds_meta.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"✓ CDS metadata saved: {output_csv}")
    print(f"  Total records: {len(cds_meta)}")
    print(f"  Human: {(cds_meta['final_label']=='human').sum()}")
    print(f"  Nonhuman: {(cds_meta['final_label']=='nonhuman').sum()}")
    
    return cds_meta

def label_protein_fasta(protein_fasta, meta_main, output_csv='protein_meta_labeled.csv'):
    """ Add labels for protein FASTA"""
    print("\n=== Processing Protein FASTA ===")
    # Prepare mapping tables.
    meta_main['id_clean'] = meta_main['id'].str.replace(r'\.\d+$', '', regex=True)
    id_to_label = dict(zip(meta_main['id_clean'], meta_main['final_label']))
    id_to_species_std = dict(zip(meta_main['id_clean'], meta_main['species_std']))
    id_to_segment = dict(zip(meta_main['id_clean'], meta_main['segment']))
    protein_records = []
    for record in SeqIO.parse(protein_fasta, 'fasta'):
        info = parse_protein_header(record)
        # Map labels via nucleotide_id.
        if info.get('nucleotide_id_clean'):
            info['species_std'] = id_to_species_std.get(info['nucleotide_id_clean'], 'Unknown')
            info['final_label'] = id_to_label.get(info['nucleotide_id_clean'], 'Unknown')
            info['segment'] = id_to_segment.get(info['nucleotide_id_clean'], 'Unknown')
        else:
            info['species_std'] = 'Unknown'
            info['final_label'] = 'Unknown'
            info['segment'] = 'Unknown'
        info['seq_length'] = len(record.seq)
        protein_records.append(info)
    
    protein_meta = pd.DataFrame(protein_records)
    protein_meta = protein_meta[protein_meta['final_label'] != 'Unknown']
    protein_meta.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"✓ Protein metadata saved: {output_csv}")
    print(f"  Total records: {len(protein_meta)}")
    print(f"  Human: {(protein_meta['final_label']=='human').sum()}")
    print(f"  Nonhuman: {(protein_meta['final_label']=='nonhuman').sum()}")

    return protein_meta

# =====================================================
# Phase 3: Feature engineering + machine learning.
# =====================================================

# --- k-mer feature extraction ---
def nucleotide_freq(seq):
    """Compute nucleotide frequencies."""
    seq = seq.upper()
    total = len(seq)
    counts = Counter(seq)
    freq = {nt: counts.get(nt, 0) / total for nt in "ACGT"}
    return freq

def bias_corrected_kmer(seq, k=2):
    """Compute bias-corrected k-mer frequencies."""
    seq = seq.upper()
    total = len(seq) - k + 1
    if total <= 0:
        return {}
    
    counts = Counter([seq[i:i+k] for i in range(total)])
    base_freq = nucleotide_freq(seq)
    
    bias_corr = {}
    for kmer, obs_count in counts.items():
        if not set(kmer) <= set("ACGT"):
            continue
        expected = 1
        for c in kmer:
            expected *= base_freq.get(c, 0)
        obs_freq = obs_count / total if total > 0 else 0
        bias_corr[kmer] = obs_freq / expected if expected > 0 else 0
    
    return bias_corr

def extract_kmer_features(nucleotide_fasta, meta_main, k_values=[2,3,4,5,6], output_csv='kmer_features.csv'):
    """Extract k-mer features."""
    print(f"\n=== Extracting k-mer features (k={k_values}) ===")
    
    meta_for_kmer = meta_main[meta_main['completeness'].isin(['complete_sequence', 'complete_cds'])].copy()
    id_to_meta = meta_for_kmer.set_index('id').to_dict('index')
    print(f"[INFO] K-mer features extracted for {len(id_to_meta)} sequences")

    features_list = []
    for record in SeqIO.parse(nucleotide_fasta, 'fasta'):
        if record.id not in id_to_meta:
            continue
        
        seq = str(record.seq)
        row = {'id': record.id}
        
        # Extract k-mer features for each k.
        for k in k_values:
            kmer_feats = bias_corrected_kmer(seq, k)
            for kmer, value in kmer_feats.items():
                row[f'{k}mer_{kmer}'] = value
        
        features_list.append(row)
    
    kmer_df = pd.DataFrame(features_list).fillna(0)
    
    # Merge with main table.
    kmer_df = kmer_df.merge(meta_main[['id', 'species_std', 'segment', 'final_label']], on='id')
    
    kmer_df.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"✓ k-mer features saved: {output_csv}")
    print(f"  Shape: {kmer_df.shape}")
    
    return kmer_df

# --- RSCU/CAI feature extraction ---
def calc_rscu(seq):
    """Compute RSCU values."""
    seq = seq.upper().replace('U','T')
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3])==3]
    codon_list = [c for c in codons if set(c) <= set("ATGC") and len(c)==3]
    obs = Counter(codon_list)
    
    aa_map = {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Y': ['TAT', 'TAC'],
        'H': ['CAT', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'C': ['TGT', 'TGC'],
        'W': ['TGG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'Stop': ['TAA', 'TAG', 'TGA'],
    }
    
    codon_table = [c for aa in aa_map.values() for c in aa]
    rscu = {}
    
    for aa, codons_for_aa in aa_map.items():
        n = sum([obs[codon] for codon in codons_for_aa])
        if n == 0:
            continue
        for codon in codons_for_aa:
            rscu[codon] = obs[codon] * len(codons_for_aa) / n if n > 0 else 0
    
    for codon in codon_table:
        if codon not in rscu:
            rscu[codon] = 0.0
    
    return rscu

class ImprovedCAI:
    def __init__(self, reference_sequences):
        self.codon_table = self._build_codon_table(reference_sequences)
        self.w_values = self._calculate_w_values()
    
    def _build_codon_table(self, reference_sequences):
        codon_count = Counter()
        for seq in reference_sequences:
            seq = seq.upper().replace('U', 'T')
            seq = ''.join([b for b in seq if b in 'ATGC'])
            if len(seq) % 3 != 0 or len(seq) < 21:
                continue
            codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
            codons = [c for c in codons if c not in ['TAA', 'TAG', 'TGA']]
            if len(codons) > 0:
                codon_count.update(codons)
        return codon_count
    
    def _calculate_w_values(self):
        genetic_code = {
            'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
            'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']
        }
        
        w_values = {}
        for aa, codons in genetic_code.items():
            if len(codons) == 1:
                w_values[codons[0]] = 1.0
            else:
                codon_freqs = {codon: self.codon_table.get(codon, 0) for codon in codons}
                max_freq = max(codon_freqs.values())
                if max_freq > 0:
                    for codon in codons:
                        w_values[codon] = codon_freqs[codon] / max_freq
                else:
                    for codon in codons:
                        w_values[codon] = 1.0 / len(codons)
        
        return w_values
    
    def calculate_cai(self, sequence):
        try:
            seq = sequence.upper().replace('U', 'T')
            seq = ''.join([b for b in seq if b in 'ATGC'])
            if len(seq) % 3 != 0 or len(seq) < 3:
                return 0.0
            
            codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
            codons = [c for c in codons if c not in ['TAA', 'TAG', 'TGA'] and set(c) <= set('ATGC')]
            
            if len(codons) == 0:
                return 0.0
            
            log_sum = 0
            valid_codons = 0
            for codon in codons:
                if codon in self.w_values and self.w_values[codon] > 0:
                    log_sum += np.log(self.w_values[codon])
                    valid_codons += 1
            
            if valid_codons == 0:
                return 0.0
            
            cai = np.exp(log_sum / valid_codons)
            return cai
        except Exception:
            return 0.0

def extract_rscu_cai_features(cds_fasta, meta_main, reference_fasta='human_HK_CDS.cleaned.fasta', 
                             output_csv='rscu_cai_features.csv', min_partial_length=600):
    """Extract RSCU/CAI from raw CDS FASTA — keep all complete_cds; 
    retain partial_cds only if they meet the length threshold."""
    print(f"\n=== Extracting RSCU/CAI from RAW FASTA (partial threshold: >{min_partial_length} nt) ===")
    #  Check CDS FASTA exists.
    if not os.path.exists(cds_fasta):
        raise FileNotFoundError(f"CDS FASTA file not found: {cds_fasta}")
    # Prepare CAI calculator.
    ref_sequences = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        seq = str(record.seq).upper().replace('U','T')
        seq = ''.join([b for b in seq if b in 'ATGC'])
        if len(seq) >= 21 and len(seq) % 3 == 0:
            ref_sequences.append(seq)
    cai_calculator = ImprovedCAI(ref_sequences)
    #  Build meta_main clean-ID mapping.
    meta_main['id_clean'] = meta_main['id'].str.replace(r'\.\d+$', '', regex=True)
    id_to_meta = meta_main.set_index('id_clean').to_dict('index')
    
    print(f"Meta main records: {len(meta_main)}")
    print(f"Unique clean IDs in meta: {len(id_to_meta)}")
    
    features_list = []
    total_records = 0
    retained_records = 0
    skip_reasons = {'no_nuc_id': 0, 'not_in_meta': 0, 'filtered_out': 0}

    #  Check CDS file has records.
    cds_record_count = 0
    for record in SeqIO.parse(cds_fasta, 'fasta'):
        cds_record_count += 1
        if cds_record_count <= 5:  # Print first 5 records as examples.
            print(f"Sample CDS record {cds_record_count}: {record.id}")
        break
    
    if cds_record_count == 0:
        print("ERROR: No records found in CDS FASTA file!")
        return pd.DataFrame()
    
    print(f"Total CDS records in file: {cds_record_count}")
    
    # Restart processing.
    for record in SeqIO.parse(cds_fasta, 'fasta'):
        total_records += 1
        if total_records <= 5:  # Verbose debug for the first 5 records.
            print(f"\nProcessing record {total_records}: {record.id}")
        
        header_info = parse_cds_header(record)
        nuc_id_clean = header_info.get('nucleotide_id_clean')
        completeness = header_info.get('completeness', 'unknown').lower()
        seq_length = len(record.seq)
        
        if total_records <= 5:
            print(f"  Header info: {header_info}")
            print(f"  Nucleotide ID clean: {nuc_id_clean}")
            print(f"  Completeness: {completeness}")
            print(f"  Sequence length: {seq_length}")
        
        if not nuc_id_clean:
            skip_reasons['no_nuc_id'] += 1
            if total_records <= 5:
                print(f"  SKIPPED: No nucleotide ID")
            continue
            
        if nuc_id_clean not in id_to_meta:
            skip_reasons['not_in_meta'] += 1
            if total_records <= 5:
                print(f"  SKIPPED: Not in meta (ID: {nuc_id_clean})")
                print(f"  Available meta IDs sample: {list(id_to_meta.keys())[:5]}")
            continue
        
        # Inclusion policy.
        if completeness == 'complete_cds':
            pass  # Keep.
        elif completeness == 'partial_cds' and seq_length >= min_partial_length:
            pass  # Keep partial_cds meeting length.
        elif completeness == 'unknown' and seq_length >= min_partial_length:
            pass  # For unknown, filter by length.
        else:
            skip_reasons['filtered_out'] += 1
            if total_records <= 5:
                print(f"  SKIPPED: Filtered out (completeness: {completeness}, length: {seq_length})")
            continue
        
        meta = id_to_meta[nuc_id_clean]
        retained_records += 1
        
        if total_records <= 5:
            print(f"  RETAINED: Adding to features list")
        
        row = {
            'cds_id': record.id,
            'nucleotide_id': meta['id'],
            'species_std': meta['species_std'],
            'segment': meta['segment'],
            'final_label': meta['final_label']
        }
        # RSCU features.
        rscu = calc_rscu(str(record.seq))
        for codon, val in rscu.items():
            row[f'RSCU_{codon}'] = val
        # CAI feature.
        row['CAI'] = cai_calculator.calculate_cai(str(record.seq))
        features_list.append(row)
    
    print(f"\nProcessing summary:")
    print(f"  Total CDS records processed: {total_records}")
    print(f"  Records retained: {retained_records}")
    print(f"  Skip reasons:")
    for reason, count in skip_reasons.items():
        print(f"    {reason}: {count}")
    # Check if features_list is empty.
    if not features_list:
        print("ERROR: No features extracted! All records were filtered out.")
        print("Check your CDS file and meta_main matching.")
        return pd.DataFrame()
    rscu_df = pd.DataFrame(features_list).fillna(0)

    # Add clean ID column
    rscu_df['nucleotide_id_clean'] = rscu_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
    
    print("rscu_df columns:", rscu_df.columns.tolist())
    print("rscu_df shape:", rscu_df.shape)
    print("rscu_df head:")
    print(rscu_df.head())
    
    rscu_df.to_csv(output_csv, index=False, encoding='utf-8')
    
    print(f"✓ RSCU/CAI features saved: {output_csv}")
    print(f"  Shape: {rscu_df.shape}")
    print(f"  Total CDS parsed: {total_records}")
    print(f"  Records matched and retained: {retained_records}")
    
    return rscu_df

# --- Protein feature extraction ---
def extract_protein_features(protein_fasta, protein_meta, meta_main, output_csv='protein_features.csv'):
    """Extract protein features via iFeatureOmegaCLI."""
    print("\n=== Extract Protein Features ===")
    if not HAS_IFEATURE:
        print("✗ iFeatureOmegaCLI not available — skipping protein feature extraction.")
        return pd.DataFrame()
    # Keep proteins derived from complete sequences.
    valid_ids = meta_main[meta_main['completeness'].isin(['complete_sequence', 'complete_cds'])]['id'].str.replace(r'\.\d+$', '', regex=True)
    protein_meta['nucleotide_id_clean'] = protein_meta['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
    protein_meta = protein_meta[protein_meta['nucleotide_id_clean'].isin(valid_ids)]
    print(f"[INFO] Protein entries retained after filtering by completeness: {len(protein_meta)}")
    # Prepare temporary FASTA.
    temp_fasta = 'temp_protein_for_features.fasta'
    protein_id_to_meta = protein_meta.set_index('protein_id').to_dict('index')
    # Write temporary FASTA (with label info).
    with open(temp_fasta, 'w', encoding='utf-8') as f:
        for record in SeqIO.parse(protein_fasta, 'fasta'):
            if record.id in protein_id_to_meta:
                meta = protein_id_to_meta[record.id]
                # Simplify header for later parsing.
                f.write(f">{record.id}|{meta['species_std']}|{meta['segment']}|{meta['final_label']}\n")
                f.write(f"{str(record.seq)}\n")
    
    # Extract descriptors.
    descriptors = ["CTDC", "CTDT", "CTDD", "CTriad", "DistancePair", "PAAC"]
    feature_dfs = {}
    
    for desc in descriptors:
        try:
            print(f"  Extracting {desc}...")
            prot = iFeatureOmegaCLI.iProtein(temp_fasta)
            prot.get_descriptor(desc)
            
            if prot.encodings is not None and len(prot.encodings.columns) > 0:
                df_temp = prot.encodings.copy()
                # Clean index; keep protein_id only.
                df_temp.index = df_temp.index.astype(str).str.extract(r'^([^|>]+)')[0].str.strip()
                df_temp = df_temp[~df_temp.index.duplicated(keep='first')]
                feature_dfs[desc] = df_temp
                print(f"    ✓ {desc}: {df_temp.shape[1]} features")
        except Exception as e:
            print(f"    ✗ {desc} failed: {e}")
    
    # Merge all descriptor sets.
    if feature_dfs:
        # Find common protein IDs.
        common_ids = sorted(set.intersection(*(set(df.index) for df in feature_dfs.values())))
        merged_df = pd.DataFrame(index=common_ids)
        
        for desc, df in feature_dfs.items():
            merged_df = pd.concat([merged_df, df.loc[common_ids]], axis=1)
        # Add metadata.
        merged_df.reset_index(inplace=True)
        merged_df.rename(columns={'index': 'protein_id'}, inplace=True)
        # Merge with protein_meta.
        protein_features = merged_df.merge(
            protein_meta[['protein_id', 'nucleotide_id', 'species_std', 'segment', 'final_label']], 
            on='protein_id'
        )
        
        protein_features.to_csv(output_csv, index=False, encoding='utf-8')
        print(f"✓ Protein features saved: {output_csv}")
        print(f"  Shape: {protein_features.shape}")
        
        # Remove temporary file.
        os.remove(temp_fasta)
        
        return protein_features
    else:
        print("✗ No protein features extracted")
        return pd.DataFrame()

# --- Feature-merging strategies ---
def merge_features(meta_main, kmer_df=None, rscu_df=None, protein_df=None):
    """Merge feature types using two strategies."""
    print("\n=== Merging features ===")
    results = {}
    
    # Strategy 1: Separate models (maximise sample size).
    if kmer_df is not None:
        results['kmer_only'] = kmer_df.copy()
        print(f"  kmer_only: {kmer_df.shape}")
    
    if rscu_df is not None:
        # For RSCU, map back to meta_main via nucleotide_id.
        rscu_with_meta = rscu_df.merge(
            meta_main[['id', 'species_std', 'segment', 'final_label']], 
            left_on='nucleotide_id', 
            right_on='id',
            suffixes=('_rscu', '')
        )
        results['rscu_only'] = rscu_with_meta
        print(f"  rscu_only: {rscu_with_meta.shape}")
    
    if protein_df is not None:
        # Proteins also map by nucleotide_id.
        # Prepare version-stripped IDs for merging.
        protein_df['nucleotide_id_clean'] = protein_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
        meta_main_subset = meta_main[['id', 'species_std', 'segment', 'final_label']].copy()
        meta_main_subset['id_clean'] = meta_main_subset['id'].str.replace(r'\.\d+$', '', regex=True)
        
        # Merge using clean_id.
        protein_with_meta = protein_df.merge(
            meta_main_subset,
            left_on='nucleotide_id_clean', 
            right_on='id_clean',
            suffixes=('_protein', '')
        )
        results['protein_only'] = protein_with_meta
        print(f"  protein_only: {protein_with_meta.shape}")

    # Strategy 2: Adaptive Intersection set (fair comparison).
    # ---------------------------------------------------------------
    # adaptively choose based on number of available sets:
    #  - 3 3 sets: intersection of all three;
    #  - 2 sets: intersection of the two;
    #  - 1 set: directly use that set;
    #  - if intersection is empty: fallback order is RSCU → k-mer → protein.

    def _ensure_clean_cols():
        # Prepare nucleotide_id_clean for RSCU/Protein
        if rscu_df is not None and 'nucleotide_id_clean' not in rscu_df.columns and 'nucleotide_id' in rscu_df.columns:
            rscu_df['nucleotide_id_clean'] = rscu_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
        if protein_df is not None and 'nucleotide_id_clean' not in protein_df.columns and 'nucleotide_id' in protein_df.columns:
            protein_df['nucleotide_id_clean'] = protein_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)

    _ensure_clean_cols()

    # Assemble available sets and their clean-ID sets
    available = {}
    if kmer_df is not None:
        # Sample IDs for k-mer come from meta_main mapping (kmer_df['id'] may include version numbers)
        kmer_ids_clean = set(
            meta_main.loc[meta_main['id'].isin(kmer_df['id']), 'id_clean']
        )
        available['kmer'] = {'df': kmer_df, 'ids': kmer_ids_clean}
    if rscu_df is not None:
        rscu_ids_clean = set(rscu_df['nucleotide_id_clean']) if 'nucleotide_id_clean' in rscu_df.columns else set()
        available['rscu'] = {'df': rscu_df, 'ids': rscu_ids_clean}
    if protein_df is not None:
        protein_ids_clean = set(protein_df['nucleotide_id_clean']) if 'nucleotide_id_clean' in protein_df.columns else set()
        available['protein'] = {'df': protein_df, 'ids': protein_ids_clean}

    # If no available set exists, return
    if not available:
        return results

    # Compute intersection of available sets (only meaningful when ≥2 sets)
    from functools import reduce
    key_order_for_base = ['kmer', 'rscu', 'protein']  # When merging, prioritize k-mer as the base (if present)
    keys = list(available.keys())

    common_ids = None
    if len(keys) >= 2:
        common_ids = reduce(lambda a, b: a & b, (available[k]['ids'] for k in keys))
    elif len(keys) == 1:
        # If only one available set, no intersection needed
        common_ids = available[keys[0]]['ids']

    def _kmer_filter_to_ids(df, id_set):
        if id_set is None:
            return df.copy()
        tmp = df.copy()
        tmp['id_clean'] = tmp['id'].str.replace(r'\.\d+$', '', regex=True)
        return tmp[tmp['id_clean'].isin(id_set)]

    def _pick_base_key():
        for k in key_order_for_base:
            if k in available:
                return k
        return keys[0]

    intersect_df = None
    if common_ids and len(common_ids) > 0:
        base_key = _pick_base_key()
        base_df = available[base_key]['df']

        
        if base_key == 'kmer':
            intersect_df = _kmer_filter_to_ids(base_df, common_ids)
        elif base_key == 'rscu':          
            intersect_df = base_df[base_df['nucleotide_id_clean'].isin(common_ids)].copy()
        else:  # protein
            intersect_df = base_df[base_df['nucleotide_id_clean'].isin(common_ids)].copy()

        # Merge the remaining features
        if 'id_clean' not in intersect_df.columns:
            if 'id' in intersect_df.columns:
                intersect_df['id_clean'] = intersect_df['id'].str.replace(r'\.\d+$', '', regex=True)
            elif 'nucleotide_id_clean' in intersect_df.columns:
                intersect_df['id_clean'] = intersect_df['nucleotide_id_clean']
        # Merge RSCU
        if 'rscu' in available and base_key != 'rscu':
            rscu_common = rscu_df[rscu_df['nucleotide_id_clean'].isin(common_ids)].copy()
            rscu_features = rscu_common.drop(columns=['species_std', 'segment', 'final_label'], errors='ignore')
            rscu_features = rscu_features.drop_duplicates(subset=['nucleotide_id_clean'])
            intersect_df = intersect_df.merge(
                rscu_features, left_on='id_clean', right_on='nucleotide_id_clean', how='left', suffixes=('', '_rscu')
            )
        # Merge protein
        if 'protein' in available and base_key != 'protein':
            protein_common = protein_df[protein_df['nucleotide_id_clean'].isin(common_ids)].copy()
            protein_features = protein_common.drop(columns=['species_std', 'segment', 'final_label'], errors='ignore')
            protein_features = protein_features.drop_duplicates(subset=['nucleotide_id_clean'])
            intersect_df = intersect_df.merge(
                protein_features, left_on='id_clean', right_on='nucleotide_id_clean', how='left', suffixes=('', '_protein')
            )
    else:
        # If intersection is empty or only one set remains: fallback to a single set (priority RSCU → k-mer → protein)
        for key in ['rscu', 'kmer', 'protein']:
            if key in available and len(available[key]['ids']) > 0:
                if key == 'kmer':
                    intersect_df = _kmer_filter_to_ids(available[key]['df'], None)
                else:
                    intersect_df = available[key]['df'].copy()
                print(f"[Fallback] 'intersect' unavailable → using {key}-only feature set as the working dataset.")
                break

    if intersect_df is not None and not intersect_df.empty:
        results['intersect'] = intersect_df
        print(f"  intersect: {intersect_df.shape}")
        if common_ids is not None:
            print(f"  Common samples: {len(common_ids)}")
    
    return results

feature_sets = {
    '3mer': lambda df: [c for c in df.columns if c.startswith('3mer_')],
    '6mer': lambda df: [c for c in df.columns if c.startswith('6mer_')],
    'rscu': lambda df: [c for c in df.columns if c.startswith('RSCU_')],
    'CTriad': lambda df: [c for c in df.columns if c.startswith('CTriad')],
    'all_protein': lambda df: [c for c in df.columns if c.startswith(('CTDC', 'CTDT', 'CTDD', 'CTriad', 'DistancePair', 'PseAAC'))],
    'mixed_3mer_rscu': lambda df: [c for c in df.columns if c.startswith(('3mer_', 'RSCU_'))],
    'mixed_3mer_rscu_ctriad': lambda df: [c for c in df.columns if c.startswith(('3mer_', 'RSCU_', 'CTriad'))]
}

individual_protein_features = {
    'CTDC': lambda df: [c for c in df.columns if c.startswith('CTDC')],
    'CTDT': lambda df: [c for c in df.columns if c.startswith('CTDT')],
    'CTDD': lambda df: [c for c in df.columns if c.startswith('CTDD')],
    'DistancePair': lambda df: [c for c in df.columns if c.startswith('DistancePair')],
    'PseAAC': lambda df: [c for c in df.columns if c.startswith('PseAAC')]
}

def get_optimized_rf_params_grouped(X, y, groups, class_weight=None, ultra_quick_mode=True):
    """
    Optimise Random Forest params 
    with group-aware 5-fold CV 
    to prevent species leakage.
    
    Args:
        X — feature matrix; 
        y — labels; 
        groups — species grouping; 
        class_weight — class weights; 
        ultra_quick_mode — minimal grid.
    
    Returns:
        dict of best params.
    """
    from sklearn.model_selection import GridSearchCV
    from sklearn.ensemble import RandomForestClassifier
    import numpy as np
    
    # Define parameter grid.
    if ultra_quick_mode:
        # Ultra-quick mode: fix n_estimators=300; tune only max_depth.
        param_grid = {
            'max_depth': [10, 20, None],
            'max_features': ['sqrt', 'log2']
        }
        # 3×2=6 combos; 5-fold → 30 fits.
    else:
        # Standard mode: broader grid.
        param_grid = {
            'n_estimators': [100, 200],
            'max_depth': [10, 20, None],
            'min_samples_split': [2, 5],
            'min_samples_leaf': [1, 2],
            'max_features': ['sqrt']   # Limit to 'sqrt' to reduce complexity
        }
    
    # Base RF model—single-threaded to avoid nested parallelism.
    rf_base = RandomForestClassifier(
        n_estimators=300,  # Use the same n_estimators for tuning and final training.
        class_weight=class_weight,
        random_state=42,
        n_jobs=1  # Avoid nested parallelism.
    )
    
    # Group-aware 5-fold (StratifiedGroupKFold if available; else GroupKFold).
    try:
        from sklearn.model_selection import StratifiedGroupKFold
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    except Exception:
        from sklearn.model_selection import GroupKFold
        cv = GroupKFold(n_splits=5) 

    # Quick check: no group overlap in the first fold.
    for tr, te in cv.split(X, y, groups):
        assert set(np.asarray(groups)[tr]).isdisjoint(set(np.asarray(groups)[te])), "Group leakage in CV"
        break  # Check only the first fold to save cost.

    # Grid search—use parallel jobs.
    gs = GridSearchCV(
        rf_base, param_grid, cv=cv, scoring='f1_macro', n_jobs=-1, verbose=0
    )
    try:
        gs.fit(X, y, groups=groups)
        best = gs.best_params_
        full_params = {
            'n_estimators': 300,
            'max_depth': best.get('max_depth', None),
            'min_samples_split': best.get('min_samples_split', 2),
            'min_samples_leaf': best.get('min_samples_leaf', 1),
            'max_features': best.get('max_features', 'sqrt')
        }
        print(f"  Best params (group-aware CV): max_depth={full_params['max_depth']}, "
              f"max_features={full_params['max_features']}, CV={gs.best_score_:.3f}")
        return full_params, float(gs.best_score_)
    except Exception as e:
        print(f"  Parameter optimization failed: {e}, using defaults")
        return {
            'n_estimators': 300, 'max_depth': None,
            'min_samples_split': 2, 'min_samples_leaf': 1, 'max_features': 'sqrt'
        }, None

def _param_file(outdir, feature_name, segment, class_weight):
    tag = "balanced" if class_weight == "balanced" else str(class_weight)
    safe = tag.replace(" ", "").replace(":", "_")
    return os.path.join(outdir, f"best_params__{feature_name}__{segment}__{safe}.json")

def save_best_params(outdir, feature_name, segment, class_weight, best_params, best_score=None):
    os.makedirs(outdir, exist_ok=True)
    payload = {"feature_set": feature_name, "segment": segment,
               "class_weight": class_weight, "best_params": best_params}
    if best_score is not None:
        payload["best_cv_score_f1_macro"] = float(best_score)
    with open(_param_file(outdir, feature_name, segment, class_weight), "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

def save_best_params_index(index_path, feature_name, segment, class_weight, best_params, best_score=None):
    os.makedirs(os.path.dirname(index_path), exist_ok=True)
    if os.path.exists(index_path):
        with open(index_path, "r", encoding="utf-8") as f:
            idx = json.load(f)
    else:
        idx = {}
    key = f"{feature_name}__{segment}__{class_weight}"
    idx[key] = {
        "feature_set": feature_name,
        "segment": segment,
        "class_weight": class_weight,
        "best_params": best_params,
        "best_cv_score_f1_macro": float(best_score) if best_score is not None else None
    }
    with open(index_path, "w", encoding="utf-8") as f:
        json.dump(idx, f, ensure_ascii=False, indent=2)

def train_test_evaluation_with_auc(df, feature_cols, segment='all', class_weight=None, optimize_params=False):
    """
    Train-test evaluation with AUC;
    parameter tuning uses stratified group-aware CV.
    """
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, recall_score, f1_score, roc_auc_score
    import numpy as np
    
    # Segment filter.
    if segment != 'all':
        df = df[df['segment'] == segment].copy()
    if len(df) < 10:
        return None
    # Prepare data.
    X = df[feature_cols].fillna(0).values
    y_labels = df['final_label'].values
    y = (y_labels == 'human').astype(int)
    # Normalise class_weight mapping.
    processed_class_weight = None
    if class_weight is not None and isinstance(class_weight, dict):
        processed_class_weight = {}
        for key, value in class_weight.items():
            if key == 'human':
                processed_class_weight[1] = value
            elif key == 'nonhuman':
                processed_class_weight[0] = value
            else:
                processed_class_weight[key] = value
    else:
        processed_class_weight = class_weight
    # Check class balance for stratification.
    unique_labels = np.unique(y)
    if len(unique_labels) < 2:
        print(f"Warning: Only one class present in segment {segment}")
        return None
    # Split data with fixed seed for reproducibility.
    try:
        X_train, X_test, y_train, y_test, species_train, species_test = train_test_split(
            X, y, df['species_std'].values,
            test_size=0.2, random_state=42, stratify=y
        )
    except ValueError:
        X_train, X_test, y_train, y_test, species_train, species_test = train_test_split(
            X, y, df['species_std'].values,
            test_size=0.2, random_state=42
        )
    # Get best params (group-aware CV).
    if optimize_params and len(X_train) >= 50:
        best_params, _ = get_optimized_rf_params_grouped(  # ← unpack: (params, cv_score).
            X_train, y_train, species_train, processed_class_weight, ultra_quick_mode=True
        )
    else:
        best_params = {
            'n_estimators': 300,
            'max_depth': None,
            'min_samples_split': 2,
            'min_samples_leaf': 1,
            'max_features': 'sqrt'  # ← add this line to align with CV params.
        }

    # Train model—use all cores.
    clf = RandomForestClassifier(
        class_weight=processed_class_weight,
        random_state=42,
        n_jobs=-1,  # Allow full parallelism during training.
        **best_params
    )
    clf.fit(X_train, y_train)
    
    # Predict.
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)
    
    # Get probability of the positive class (human=1).
    if len(clf.classes_) == 2 and 1 in clf.classes_:
        human_idx = list(clf.classes_).index(1)
        y_proba_human = y_proba[:, human_idx]
    else:
        y_proba_human = np.full(len(y_test), np.nan)
    
    # Compute metrics.
    accuracy = accuracy_score(y_test, y_pred)
    
    # Compute recall/F1.
    try:
        recall_human = recall_score(y_test, y_pred, pos_label=1, zero_division=0)
        recall_nonhuman = recall_score(y_test, y_pred, pos_label=0, zero_division=0)
        f1_human = f1_score(y_test, y_pred, pos_label=1, zero_division=0)
        f1_nonhuman = f1_score(y_test, y_pred, pos_label=0, zero_division=0)
    except Exception as e:
        print(f"Warning: Error calculating metrics: {e}")
        recall_human = recall_nonhuman = f1_human = f1_nonhuman = np.nan
    
    # Compute AUC.
    try:
        if len(np.unique(y_test)) > 1 and not np.any(np.isnan(y_proba_human)):
            auc = roc_auc_score(y_test, y_proba_human)
        else:
            auc = np.nan
    except Exception as e:
        print(f"Warning: AUC calculation failed: {e}")
        auc = np.nan
    
    return {
        'accuracy': accuracy,
        'recall_human': recall_human,
        'recall_nonhuman': recall_nonhuman,
        'f1_human': f1_human,
        'f1_nonhuman': f1_nonhuman,
        'auc': auc,
        'n_samples': len(df),
        'n_features': len(feature_cols),
        'optimized_params': best_params if optimize_params else None
    }

def leave_one_virus_out_cv_with_segments(df, feature_cols, class_weight=None, segments=['all'],
    max_species=None, optimize_params=False, feature_name=None, outdir="results"):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, recall_score, f1_score, roc_auc_score
    import numpy as np
    import pandas as pd 
    
    # Normalise class_weight mapping
    processed_class_weight = None
    if class_weight is not None and isinstance(class_weight, dict):
        processed_class_weight = {}
        for key, value in class_weight.items():
            if key == 'human':
                processed_class_weight[1] = value
            elif key == 'nonhuman':
                processed_class_weight[0] = value
            else:
                processed_class_weight[key] = value
    else:
        processed_class_weight = class_weight
    
    results = []
    segment_pooled_results = {}
    
    for segment in segments:
        print(f"\n--- Processing segment: {segment} ---")
        
        # Segment filter.
        if segment != 'all':
            seg_df = df[df['segment'] == segment].copy()
        else:
            seg_df = df.copy()
            
        if len(seg_df) < 10:
            print(f"Skipping segment {segment}: insufficient data ({len(seg_df)} samples)")
            continue
            
        seg_species = seg_df['species_std'].unique()
        if max_species:
            seg_species = seg_species[:max_species]
        
        print(f"Processing {len(seg_species)} species in segment {segment}")
        
        # Optimise params at the segment level
        segment_best_params = None
        if optimize_params and len(seg_df) >= 50:
            print(f"  Optimizing parameters for segment {segment} using group-aware 5-fold CV...")
            X_seg = seg_df[feature_cols].fillna(0).values
            y_seg = (seg_df['final_label'].values == 'human').astype(int)
            groups_seg = seg_df['species_std'].values
            
            # Optimise on the whole segment with group-aware CV.
            # Ensures a species never appears in both train and validation.
            segment_best_params, seg_cv = get_optimized_rf_params_grouped(
                X_seg, y_seg, groups_seg, processed_class_weight, ultra_quick_mode=True
            )
            # Persist results (if feature_name/outdir provided).
            if feature_name is not None:
                save_best_params_index(os.path.join(outdir, "best_params_index.json"),
                       feature_name, segment, str(class_weight),
                       segment_best_params, seg_cv)
        else:
            segment_best_params = {'n_estimators': 300, 'max_depth': None,
                                'min_samples_split': 2, 'min_samples_leaf': 1, 'max_features': 'sqrt'}
            seg_cv = None
        print(f"  Using params for segment={segment}: {segment_best_params} ; CV={seg_cv}")
        
        # Init pooled containers for this segment.
        pooled_y_true = []
        pooled_y_pred = []
        pooled_y_proba = []
        pooled_test_info = []
        
        # Run LOVO with fixed segment-level params
        for test_species in seg_species:
            # Split data.
            train_df = seg_df[seg_df['species_std'] != test_species]
            test_df = seg_df[seg_df['species_std'] == test_species]
            
            # Quick check: the held-out species is absent from training.
            assert test_species not in set(train_df['species_std']), "LOVO leakage"

            if len(train_df) < 10 or len(test_df) < 1:
                continue
            # Prepare data.
            X_train = train_df[feature_cols].fillna(0).values
            y_train_labels = train_df['final_label'].values
            y_train = (y_train_labels == 'human').astype(int)
            
            X_test = test_df[feature_cols].fillna(0).values
            y_test_labels = test_df['final_label'].values
            y_test = (y_test_labels == 'human').astype(int)
            
            # Verify both classes exist in training.
            if len(np.unique(y_train)) < 2:
                continue
                    
            # Train model (using segment-level tuned params).
            clf = RandomForestClassifier(
                class_weight=processed_class_weight,
                random_state=42,
                n_jobs=-1,  # Parallelism allowed inside LOVO.
                **segment_best_params  # Use fixed tuned params.
            )
            clf.fit(X_train, y_train)
            
            # Predict.
            y_pred = clf.predict(X_test)
            y_proba = clf.predict_proba(X_test)
            
            # Get probability for the human class.
            if len(clf.classes_) == 2 and 1 in clf.classes_:
                human_idx = list(clf.classes_).index(1)
                y_proba_human = y_proba[:, human_idx]
            else:
                y_proba_human = np.full(len(y_test), np.nan)
            
            # Accumulate into pooled containers.
            pooled_y_true.extend(y_test)
            pooled_y_pred.extend(y_pred)
            pooled_y_proba.extend(y_proba_human)
            
            # Record test-sample details.
            for sample_idx, (true_label, pred_label, prob) in enumerate(zip(y_test, y_pred, y_proba_human)):
                pooled_test_info.append({
                    'segment': segment,
                    'test_species': test_species,
                    'sample_idx': sample_idx,
                    'y_true': true_label,
                    'y_pred': pred_label,
                    'y_proba_human': prob,
                    'true_label_str': 'human' if true_label == 1 else 'nonhuman',
                    'pred_label_str': 'human' if pred_label == 1 else 'nonhuman'
                })
            
            # Compute metrics.
            accuracy = accuracy_score(y_test, y_pred)
            
            # AUC calculation.
            try:
                if len(np.unique(y_test)) > 1 and not np.any(np.isnan(y_proba_human)):
                    auc = roc_auc_score(y_test, y_proba_human)
                else:
                    auc = np.nan
            except Exception:
                auc = np.nan
            
            # Per-class recall and F1-score.
            try:
                recall_per_class = recall_score(y_test, y_pred, average=None, labels=[0, 1], zero_division=0)
                f1_per_class = f1_score(y_test, y_pred, average=None, labels=[0, 1], zero_division=0)
                
                if len(recall_per_class) >= 2:
                    recall_nonhuman = float(recall_per_class[0])
                    recall_human = float(recall_per_class[1])
                else:
                    recall_nonhuman = np.nan
                    recall_human = float(recall_per_class[0]) if len(recall_per_class) > 0 else np.nan
                
                if len(f1_per_class) >= 2:
                    f1_nonhuman = float(f1_per_class[0])
                    f1_human = float(f1_per_class[1])
                else:
                    f1_nonhuman = np.nan
                    f1_human = float(f1_per_class[0]) if len(f1_per_class) > 0 else np.nan
                
                recall_macro = float(recall_score(y_test, y_pred, average='macro', zero_division=0))
                f1_macro = float(f1_score(y_test, y_pred, average='macro', zero_division=0))
                recall_micro = float(recall_score(y_test, y_pred, average='micro', zero_division=0))
                f1_micro = float(f1_score(y_test, y_pred, average='micro', zero_division=0))
                
            except Exception as e:
                print(f"Warning: Error calculating metrics for {test_species}: {e}")
                recall_nonhuman = recall_human = recall_macro = recall_micro = np.nan
                f1_nonhuman = f1_human = f1_macro = f1_micro = np.nan
            
            # Count class distribution in the test set.
            n_human = int(np.sum(y_test == 1))
            n_nonhuman = int(np.sum(y_test == 0))
            
            if n_human > n_nonhuman:
                target_type = 'human'
            elif n_nonhuman > n_human:
                target_type = 'nonhuman'
            else:
                target_type = 'mixed'
            
            result = {
                'segment': str(segment),
                'test_virus': str(test_species),
                'species': str(test_species),
                'accuracy': float(accuracy),
                'auc': float(auc) if not np.isnan(auc) else np.nan,
                'recall_human': float(recall_human) if not np.isnan(recall_human) else np.nan,
                'recall_nonhuman': float(recall_nonhuman) if not np.isnan(recall_nonhuman) else np.nan,
                'f1_human': float(f1_human) if not np.isnan(f1_human) else np.nan,
                'f1_nonhuman': float(f1_nonhuman) if not np.isnan(f1_nonhuman) else np.nan,
                'recall_macro': float(recall_macro) if not np.isnan(recall_macro) else np.nan,
                'recall_micro': float(recall_micro) if not np.isnan(recall_micro) else np.nan,
                'f1_macro': float(f1_macro) if not np.isnan(f1_macro) else np.nan,
                'f1_micro': float(f1_micro) if not np.isnan(f1_micro) else np.nan,
                'recall': float(recall_macro) if not np.isnan(recall_macro) else np.nan,
                'f1_score': float(f1_macro) if not np.isnan(f1_macro) else np.nan,
                'n_samples': int(len(y_test)),
                'n_human': int(n_human),
                'n_nonhuman': int(n_nonhuman),
                'target_type': str(target_type),
                'class_weight': str(class_weight) if class_weight is not None else 'None',
                'result_type': 'per_virus',
                'optimized_params': segment_best_params if optimize_params else None
            }
            
            results.append(result)
        
        # Compute pooled metrics for this segment.
        if pooled_y_true and len(pooled_y_true) > 0:
            pooled_y_true = np.array(pooled_y_true)
            pooled_y_pred = np.array(pooled_y_pred)
            pooled_y_proba = np.array(pooled_y_proba)
            
            segment_pooled_results[segment] = {
                'y_true': pooled_y_true,
                'y_pred': pooled_y_pred,
                'y_proba': pooled_y_proba,
                'test_info': pooled_test_info
            }
            
            print(f"Segment {segment}: {len(pooled_y_true)} total test samples")
            print(f"  Class distribution - Human: {np.sum(pooled_y_true)}, Non-human: {np.sum(pooled_y_true == 0)}")
            
            try:
                pooled_accuracy = float(accuracy_score(pooled_y_true, pooled_y_pred))
                
                valid_proba_mask = ~np.isnan(pooled_y_proba)
                if np.sum(valid_proba_mask) > 0 and len(np.unique(pooled_y_true[valid_proba_mask])) > 1:
                    pooled_auc = float(roc_auc_score(pooled_y_true[valid_proba_mask], pooled_y_proba[valid_proba_mask]))
                else:
                    pooled_auc = np.nan
                
                pooled_recall_per_class = recall_score(pooled_y_true, pooled_y_pred, average=None, labels=[0, 1], zero_division=0)
                pooled_f1_per_class = f1_score(pooled_y_true, pooled_y_pred, average=None, labels=[0, 1], zero_division=0)
                
                pooled_recall_nonhuman = float(pooled_recall_per_class[0])
                pooled_recall_human = float(pooled_recall_per_class[1])
                pooled_f1_nonhuman = float(pooled_f1_per_class[0])
                pooled_f1_human = float(pooled_f1_per_class[1])
                
                pooled_recall_macro = float(recall_score(pooled_y_true, pooled_y_pred, average='macro', zero_division=0))
                pooled_f1_macro = float(f1_score(pooled_y_true, pooled_y_pred, average='macro', zero_division=0))
                pooled_recall_micro = float(recall_score(pooled_y_true, pooled_y_pred, average='micro', zero_division=0))
                pooled_f1_micro = float(f1_score(pooled_y_true, pooled_y_pred, average='micro', zero_division=0))
                
                pooled_result = {
                    'segment': str(segment),
                    'test_virus': f'POOLED_{segment}',
                    'species': f'POOLED_{segment}',
                    'accuracy': float(pooled_accuracy),
                    'auc': float(pooled_auc) if not np.isnan(pooled_auc) else np.nan,
                    'recall_human': float(pooled_recall_human),
                    'recall_nonhuman': float(pooled_recall_nonhuman),
                    'f1_human': float(pooled_f1_human),
                    'f1_nonhuman': float(pooled_f1_nonhuman),
                    'recall_macro': float(pooled_recall_macro),
                    'recall_micro': float(pooled_recall_micro),
                    'f1_macro': float(pooled_f1_macro),
                    'f1_micro': float(pooled_f1_micro),
                    'recall': float(pooled_recall_macro),
                    'f1_score': float(pooled_f1_macro),
                    'n_samples': int(len(pooled_y_true)),
                    'n_human': int(np.sum(pooled_y_true)),
                    'n_nonhuman': int(np.sum(pooled_y_true == 0)),
                    'target_type': 'mixed',
                    'class_weight': str(class_weight) if class_weight is not None else 'None',
                    'result_type': 'pooled',
                    'optimized_params': segment_best_params if optimize_params else None
                }
                
                results.append(pooled_result)
                
                print(f"✓ Pooled results for {segment}:")
                print(f"  Accuracy: {pooled_accuracy:.3f}")
                print(f"  AUC: {pooled_auc:.3f}")
                print(f"  Recall (macro): {pooled_recall_macro:.3f}")
                print(f"  F1 (macro): {pooled_f1_macro:.3f}")
                
            except Exception as e:
                print(f"✗ Error calculating pooled metrics for {segment}: {e}")
    
    # Aggregate results.
    lovo_df = pd.DataFrame(results)
    
    if lovo_df.empty:
        print("Warning: No LOVO results generated")
        return lovo_df, pd.DataFrame(), {}
    
    print(f"\n=== Data Integrity Check ===")
    critical_columns = ['species', 'recall_micro', 'f1_micro', 'recall', 'f1_score', 'target_type']
    for col in critical_columns:
        if col in lovo_df.columns:
            null_count = lovo_df[col].isna().sum()
            total_count = len(lovo_df)
            print(f"{col}: {null_count}/{total_count} null values")
    
    per_virus_df = lovo_df[lovo_df['result_type'] == 'per_virus']
    summary_results = []
    
    for segment in segments:
        seg_results = per_virus_df[per_virus_df['segment'] == segment]
        if not seg_results.empty:
            summary_results.append({
                'segment': segment,
                'n_viruses': len(seg_results),
                'avg_accuracy': seg_results['accuracy'].mean(),
                'std_accuracy': seg_results['accuracy'].std(),
                'avg_auc': seg_results['auc'].mean(),
                'std_auc': seg_results['auc'].std(),
                'avg_recall_human': seg_results['recall_human'].mean(),
                'avg_recall_nonhuman': seg_results['recall_nonhuman'].mean(),
                'avg_f1_human': seg_results['f1_human'].mean(),
                'avg_f1_nonhuman': seg_results['f1_nonhuman'].mean(),
                'avg_recall_macro': seg_results['recall_macro'].mean(),
                'avg_f1_macro': seg_results['f1_macro'].mean(),
                'std_recall_macro': seg_results['recall_macro'].std(),
                'std_f1_macro': seg_results['f1_macro'].std(),
                'n_human_dominant': len(seg_results[seg_results['target_type'] == 'human']),
                'n_nonhuman_dominant': len(seg_results[seg_results['target_type'] == 'nonhuman']),
            })
    
    summary_df = pd.DataFrame(summary_results)
    
    pooled_df = lovo_df[lovo_df['result_type'] == 'pooled']
    print(f"\n✓ LOVO evaluation completed:")
    print(f"  Per-virus results: {len(per_virus_df)}")
    print(f"  Pooled results: {len(pooled_df)}")
    print(f"  Segments: {len(summary_df)}")
    
    return lovo_df, summary_df, segment_pooled_results

def create_short_feature_names(df, col_name='feature_set'):
    if col_name not in df.columns: return df
    mapping = create_unified_feature_mapping(df[col_name].dropna().unique())
    short_col_name = f"{col_name}_short"
    df[short_col_name] = df[col_name].map(mapping).fillna(df[col_name])
    print(f"✓ Created '{short_col_name}' column for plotting.")
    return df

def run_ml_evaluation(merged_features, output_dir='results', run_lovo=False, max_species=None, 
                     optimize_params=False):
    import os
    import pandas as pd
    os.makedirs(output_dir, exist_ok=True)
    print("\n=== Machine Learning Evaluation ===")
    feature_sets = {
        '3mer': lambda df: [col for col in df.columns if col.startswith('3mer_')],
        '6mer': lambda df: [col for col in df.columns if col.startswith('6mer_')],
        'rscu': lambda df: [col for col in df.columns if col.startswith('RSCU_')],
        'CTriad': lambda df: [col for col in df.columns if col.startswith('CTriad')],
        'all_protein': lambda df: [col for col in df.columns if any(col.startswith(p) for p in ['CTDC', 'CTDT', 'CTDD', 'CTriad', 'DistancePair', 'PseAAC'])],
        'mixed_3mer_rscu': lambda df: [col for col in df.columns if col.startswith('3mer_') or col.startswith('RSCU_')],
        'mixed_3mer_rscu_ctriad': lambda df: [col for col in df.columns if col.startswith('3mer_') or col.startswith('RSCU_') or col.startswith('CTriad')]
    }
    
    if optimize_params:
        print("✓ Parameter optimization ENABLED (group-aware 5-fold CV)")
        print("  Strategy: Optimize once per segment using StratifiedGroupKFold")
    else:
        print("✗ Parameter optimization DISABLED (using default params)")
        
    # Use the intersect dataset.
    if 'intersect' not in merged_features or merged_features['intersect'] is None:
        print("Error: 'intersect' dataset not found. Cannot proceed with evaluation.")
        return None, None
        
    intersect_df = merged_features['intersect']
    print(f"Using intersect dataset with {len(intersect_df)} samples")

    # Set class weights.
    class_weights = [
        None,
        'balanced',
        {'human': 5, 'nonhuman': 1}
    ]

    segments = ['all', 'S', 'M', 'L']
        
    # 1. Train-test evaluation.
    print("\n--- Train-Test Evaluation ---")
    train_test_results = []
        
    for feat_name, feat_func in feature_sets.items():
        feat_cols = feat_func(intersect_df)
        if len(feat_cols) < 5:
            print(f"[Skip] {feat_name}: too few features ({len(feat_cols)})")
            continue
            
        print(f"\nEvaluating {feat_name} ({len(feat_cols)} features)...")
            
        for segment in segments:
            for cw in class_weights:
                try:
                    result = train_test_evaluation_with_auc(
                        intersect_df, feat_cols, segment, cw, 
                        optimize_params=optimize_params
                    )
                    
                    if result:
                        result.update({
                            'feature_set': feat_name,
                            'segment': segment,
                            'class_weight': str(cw)
                        })
                        train_test_results.append(result)
                except Exception as e:
                    print(f"[Error] Train-test failed for {feat_name} ({segment}, {cw}): {e}")
    
    # Save train–test results.
    train_test_df = pd.DataFrame(train_test_results)
    train_test_df.to_csv(os.path.join(output_dir, 'train_test_results.csv'), 
                        encoding='utf-8', index=False)
    print(f"\n✓ Train-test results saved to {output_dir}/train_test_results.csv")
        
    # 2. LOVO evaluation.
    print("\n--- LOVO Cross-Validation with Pooled Results ---")
    if run_lovo:
        lovo_results = []
        all_pooled_data = {}
        
        for feat_name, feat_func in feature_sets.items():
            feat_cols = feat_func(intersect_df)
            if len(feat_cols) < 5:
                continue
                
            print(f"\nLOVO evaluation for {feat_name}...")
            feat_pooled_data = {}
            
            for cw in class_weights:
                try:
                    lovo_df, sample_df, pooled_data = leave_one_virus_out_cv_with_segments(
                        intersect_df, feat_cols, cw, segments,
                        max_species=max_species, optimize_params=optimize_params,
                        feature_name=feat_name, outdir=output_dir
                    )

                    if not lovo_df.empty:
                        lovo_df['feature_set'] = feat_name
                        lovo_df['class_weight'] = str(cw)
                        lovo_results.append(lovo_df)
                        # Store aggregated data.
                        feat_pooled_data[str(cw)] = pooled_data
                        
                except Exception as e:
                    print(f"[Error] LOVO failed for {feat_name} ({cw}): {e}")
            all_pooled_data[feat_name] = feat_pooled_data

        if lovo_results:
            # Combine all LOVO outputs.
            all_lovo = pd.concat(lovo_results, ignore_index=True)
            # Save full results.
            all_lovo.to_csv(os.path.join(output_dir, 'lovo_results.csv'), 
                           index=False, encoding='utf-8')
            print(f"\n✓ Complete LOVO results saved to {output_dir}/lovo_results.csv")
            # Save per-virus and pooled results separately.
            per_virus_results = all_lovo[all_lovo['result_type'] == 'per_virus']
            pooled_results = all_lovo[all_lovo['result_type'] == 'pooled']
            
            per_virus_results.to_csv(os.path.join(output_dir, 'lovo_per_virus_results.csv'), 
                                   index=False, encoding='utf-8')
            pooled_results.to_csv(os.path.join(output_dir, 'lovo_pooled_results.csv'), 
                                 index=False, encoding='utf-8')
            
            print(f"✓ Per-virus results: {len(per_virus_results)} records")
            print(f"✓ Pooled results: {len(pooled_results)} records")
            # Optionally save per-sample test details.
            if all_pooled_data:
                detailed_samples = []
                for feat_name, feat_data in all_pooled_data.items():
                    for cw, cw_data in feat_data.items():
                        for segment, seg_data in cw_data.items():
                            for sample_info in seg_data['test_info']:
                                sample_info.update({
                                    'feature_set': feat_name,
                                    'class_weight': cw
                                })
                                detailed_samples.append(sample_info)
                
                if detailed_samples:
                    detailed_df = pd.DataFrame(detailed_samples)
                    detailed_df.to_csv(os.path.join(output_dir, 'lovo_detailed_samples.csv'), 
                                     index=False, encoding='utf-8')
                    print(f"✓ Detailed sample results: {len(detailed_df)} records")
            # Generate comparison report.
            print(f"\n--- Performance Comparison Report ---")
            # Show best configurations.
            if optimize_params:
                print("\nOptimized Parameters Used:")
                param_cols = [col for col in all_lovo.columns if 'param' in col.lower()]
                if 'optimized_params' in all_lovo.columns:
                    unique_params = all_lovo[all_lovo['optimized_params'].notna()]['optimized_params'].apply(str).unique()
                    for i, params in enumerate(unique_params[:3]):  # Show top 3 only.
                        print(f"  Config {i+1}: {params}")
            print("\nPer-virus vs Pooled Results:")
            
            for feat_name in all_lovo['feature_set'].unique():
                print(f"\n{feat_name.upper()}:")
                
                feat_data = all_lovo[all_lovo['feature_set'] == feat_name]
                
                for cw in feat_data['class_weight'].unique():
                    cw_data = feat_data[feat_data['class_weight'] == cw]
                    
                    per_virus_subset = cw_data[cw_data['result_type'] == 'per_virus']
                    pooled_subset = cw_data[cw_data['result_type'] == 'pooled']
                    
                    if not per_virus_subset.empty and not pooled_subset.empty:
                        pv_acc_mean = per_virus_subset['accuracy'].mean()
                        pv_f1_mean = per_virus_subset['f1_macro'].mean()
                        
                        pooled_acc_mean = pooled_subset['accuracy'].mean()
                        pooled_f1_mean = pooled_subset['f1_macro'].mean()
                        
                        print(f"  Class Weight {cw}:")
                        print(f"    Per-virus avg: Acc={pv_acc_mean:.3f}, F1={pv_f1_mean:.3f}")
                        print(f"    Pooled avg:    Acc={pooled_acc_mean:.3f}, F1={pooled_f1_mean:.3f}")
                        print(f"    Difference:    Acc={pooled_acc_mean - pv_acc_mean:+.3f}, F1={pooled_f1_mean - pv_f1_mean:+.3f}")
        else:
            all_lovo = None
            print("No LOVO results generated.")
    else:
        # If evaluation is skipped, try loading existing results.
        try:
            all_lovo = pd.read_csv(os.path.join(output_dir, 'lovo_results.csv'))
            all_lovo = create_short_feature_names(all_lovo, 'feature_set')
            print("✓ Loaded existing LOVO results.")
        except Exception as e:
            all_lovo = None
            print(f"✗ Failed to load existing LOVO results: {e}")

    return train_test_df if not train_test_df.empty else None, \
           all_lovo if all_lovo is not None and not all_lovo.empty else None

from sklearn.model_selection import StratifiedGroupKFold

def sanity_check_group_cv(X, y, groups):
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    for k, (tr, te) in enumerate(cv.split(X, y, groups)):
        g_tr, g_te = set(groups[tr]), set(groups[te])
        inter = g_tr.intersection(g_te)
        assert len(inter) == 0, f"Fold {k}: group leakage: {inter}"
    print("✓ No group leakage across folds.")


def debug_plot_environment():
    print(f"Matplotlib backend: {matplotlib.get_backend()}")
    print(f"Matplotlib version: {matplotlib.__version__}")
    print(f"Current working directory: {os.getcwd()}")
    
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([1, 2, 3], [1, 4, 2])
        ax.set_title("Test Plot")
        
        os.makedirs('results', exist_ok=True)
        test_file = 'results/test_plot.png'
        plt.savefig(test_file, dpi=200, bbox_inches='tight')
        plt.close()
        
        if os.path.exists(test_file):
            print(f"✓ Test plot saved successfully: {test_file}")
            print(f"File size: {os.path.getsize(test_file)} bytes")
            return True
        else:
            print("✗ Test plot file not created")
            return False
    except Exception as e:
        print(f"✗ Test plot failed: {e}")
        return False

def _metric_available(df, col):
    return (col in df.columns) and pd.api.types.is_numeric_dtype(df[col]) and df[col].notna().any()

def parse_weight(weight):
    if weight is None or (isinstance(weight, float) and pd.isna(weight)):
        return 'None'
    s = str(weight).strip().lower()
    if s in {'none','nan','na'}:
        return 'None'
    if 'balanced' in s:
        return 'Balanced'
    if re.search(r"human['\"]?\s*:\s*5", s) and re.search(r"nonhuman['\"]?\s*:\s*1", s):
        return '5:1'
    return s

def create_unified_feature_mapping(feature_sets):
    mapping = {}
    for fs in sorted(feature_sets):
        low = fs.lower()
        if 'mixed_3mer_rscu_ctriad' in low: mapping[fs] = 'Mix3'
        elif 'mixed_3mer_rscu' in low:      mapping[fs] = 'Mix2'
        elif 'all_protein' in low:          mapping[fs] = 'AllProt'
        elif 'ctriad' in low:               mapping[fs] = 'CTriad'
        elif '6mer' in low:                 mapping[fs] = '6mer'
        elif '3mer' in low:                 mapping[fs] = '3mer'
        elif 'rscu' in low:                 mapping[fs] = 'RSCU'
        else:                               mapping[fs] = fs[:6] if len(fs)>6 else fs
    return mapping

def plot_umap(feature_df, feature_cols, title, output_file):
    print(f"\n--- Starting UMAP: {title} ---")
    # Data validation
    if feature_df is None or feature_df.empty:
        print(f"Skip UMAP: input dataframe is empty for {title}")
        return False
    if not feature_cols or len(feature_cols) < 2:
        print(f"Skip UMAP: insufficient features for {title} (got {len(feature_cols) if feature_cols else 0} features)")
        return False
    if len(feature_df) < 10:
        print(f"Skip UMAP: insufficient samples for {title} (got {len(feature_df)} samples)")
        return False
    print(f"Input data shape: {feature_df.shape}")
    print(f"Feature columns: {len(feature_cols)}")
    # Check label columns.
    if 'final_label' not in feature_df.columns:
        print(f"✗ Missing 'final_label' column for {title}")
        return False
    unique_labels = feature_df['final_label'].unique()
    print(f"Unique labels: {unique_labels}")
    
    try:
        # Data preprocessing.
        X = feature_df[feature_cols].fillna(0)
        # Check for inf/NaN values.
        if np.isinf(X.values).any():
            print("Warning: Infinite values detected, replacing with 0")
            X = X.replace([np.inf, -np.inf], 0)
        
        if X.isnull().any().any():
            print("Warning: NaN values still present after fillna")
            X = X.fillna(0)
        # Attempt to import UMAP.
        try:
            import umap
            print("✓ UMAP imported successfully")
        except ImportError:
            print("✗ UMAP not available, skipping")
            return False
        # Run UMAP dimensionality reduction.
        print("Performing UMAP reduction...")
        reducer = umap.UMAP(random_state=42, n_neighbors=15, min_dist=0.1)
        embedding = reducer.fit_transform(X.values)
        print(f"UMAP embedding shape: {embedding.shape}")
        plt.figure(figsize=(8, 6))
        # Plot points by label.
        colors = {'human': '#e24a33', 'nonhuman': '#348abd'}
        for label in unique_labels:
            if label in colors:
                idx = (feature_df['final_label'] == label)
                count = idx.sum()
                print(f"Plotting {count} points for label '{label}'")
                
                if count > 0:
                    plt.scatter(embedding[idx, 0], embedding[idx, 1], 
                               label=f'{label} (n={count})', s=30, alpha=0.7, c=colors[label])
        
        plt.title(title, fontsize=12)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True) # Ensure directory exists.

        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        plt.close()
        
        # Verify output file created.
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            timestamp = datetime.fromtimestamp(os.path.getmtime(output_file))
            print(f"✓ UMAP saved: {output_file}")
            print(f"  File size: {file_size} bytes")
            print(f"  Timestamp: {timestamp}")
            return True
        else:
            print(f"✗ UMAP file not created: {output_file}")
            return False
            
    except Exception as e:
        print(f"✗ UMAP failed for {title}: {e}")
        import traceback
        traceback.print_exc()
        return False

def plot_performance_comparison_single(results_df, output_file_prefix):
    print(f"\n--- Starting Performance Comparison ---")
    
    if results_df.empty:
        print("[Warning] Performance results are empty, skipping comparison plot.")
        return False
    
    print(f"Results data shape: {results_df.shape}")
    print(f"Columns: {results_df.columns.tolist()}")
    print(f"Segments: {results_df['segment'].unique()}")
    print(f"Feature sets: {results_df['feature_set'].unique()}")
    print(f"Class weights: {results_df['class_weight'].unique()}")
    
    weight_labels = {
        'None': 'None',
        "{'human': 2, 'nonhuman': 1}": '2:1',
        "{'human': 3, 'nonhuman': 1}": '3:1',
        "{'human': 5, 'nonhuman': 1}": '5:1',
        'balanced': 'bal.'
    }
    
    segments = results_df['segment'].unique()
    success_count = 0
    
    for seg in segments:
        print(f"\nProcessing segment: {seg}")
        seg_data = results_df[results_df['segment'] == seg]
        if seg_data.empty:
            print(f"No data for segment {seg}")
            continue
        
        try:
            plt.figure(figsize=(12, 6))
            
            feature_sets = seg_data['feature_set'].unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(feature_sets)))
            
            plot_created = False
            
            for i, feature_set in enumerate(feature_sets):
                subset = seg_data[seg_data['feature_set'] == feature_set]
                if subset.empty:
                    continue
                
                x_values = []
                y_values = []
                
                # Process class-weight data.
                for class_weight in subset['class_weight'].unique():
                    weight_data = subset[subset['class_weight'] == class_weight]
                    if not weight_data.empty and 'accuracy' in weight_data.columns:
                        x_label = weight_labels.get(str(class_weight), str(class_weight))
                        accuracy_val = weight_data['accuracy'].iloc[0]
                        
                        if pd.notna(accuracy_val):
                            x_values.append(x_label)
                            y_values.append(accuracy_val)
                
                if x_values and y_values:
                    marker = 's' if feature_set == 'protein' else 'o'
                    linewidth = 3 if feature_set == 'protein' else 2
                    markersize = 8 if feature_set == 'protein' else 6
                    color = '#FF6B6B' if feature_set == 'protein' else colors[i]
                    
                    plt.plot(x_values, y_values, marker=marker, label=feature_set, 
                            color=color, linewidth=linewidth, markersize=markersize)
                    plot_created = True
                    print(f"  Plotted {feature_set}: {len(x_values)} points")
            
            if plot_created:
                plt.title(f"Accuracy by Class Weight (Segment: {seg})", fontsize=14)
                plt.xlabel('Class Weight', fontsize=12)
                plt.ylabel('Accuracy', fontsize=12)
                plt.xticks(rotation=0, ha='center')
                plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                
                output_file = f'{output_file_prefix}_{seg}.png'
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                
                # Verify file creation
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file)
                    timestamp = datetime.fromtimestamp(os.path.getmtime(output_file))
                    print(f"✓ Performance comparison saved: {output_file}")
                    print(f"  File size: {file_size} bytes")
                    print(f"  Timestamp: {timestamp}")
                    success_count += 1
                else:
                    print(f"✗ Performance comparison file not created: {output_file}")
            else:
                print(f"No valid data to plot for segment {seg}")
                plt.close()
                
        except Exception as e:
            print(f"✗ Error creating performance plot for segment {seg}: {e}")
            import traceback
            traceback.print_exc()
            plt.close()
    
    return success_count > 0

def plot_lovo_boxplot(lovo_df, output_file):

    print(f"\n--- Starting Rubost LOVO Boxplot ---")
    
    if lovo_df is None or lovo_df.empty:
        print("Empty LOVO data for boxplot")
        return False
    
    print(f"LOVO data shape: {lovo_df.shape}")
    
    if 'result_type' in lovo_df.columns:
        type_counts = lovo_df['result_type'].value_counts(dropna=False)
        print(f"Input result types: {type_counts.to_dict()}")
    
    try:
        lovo_data = lovo_df.copy()
        lovo_data['weight_label'] = lovo_data['class_weight'].apply(parse_weight)
        
        # Split data.
        if 'result_type' in lovo_data.columns:
            per_virus_data = lovo_data[lovo_data['result_type'] == 'per_virus'].copy()
            pooled_data = lovo_data[lovo_data['result_type'] == 'pooled'].copy()
        else:
            per_virus_data = lovo_data.copy()
            pooled_data = pd.DataFrame()
        
        print(f"✓ Data separation: Per-virus={len(per_virus_data)}, Pooled={len(pooled_data)}")
        
        if per_virus_data.empty:
            print("No per-virus data available")
            return False
        
        # Unified feature-set label mapping.
        def create_unified_feature_mapping(feature_sets):
            feature_mapping = {}
            for fs in sorted(feature_sets):
                if 'mixed_3mer_rscu_ctriad' in fs.lower():
                    feature_mapping[fs] = 'Mix3'
                elif 'mixed_3mer_rscu' in fs.lower():
                    feature_mapping[fs] = 'Mix2'  
                elif 'all_protein' in fs.lower():
                    feature_mapping[fs] = 'AllProt'
                elif 'ctriad' in fs.lower():
                    feature_mapping[fs] = 'CTriad'
                elif '3mer' == fs.lower():
                    feature_mapping[fs] = '3mer'
                elif '6mer' == fs.lower():
                    feature_mapping[fs] = '6mer'
                elif 'rscu' == fs.lower():
                    feature_mapping[fs] = 'RSCU'
                else:
                    feature_mapping[fs] = fs[:6] if len(fs) > 6 else fs
            return feature_mapping
        
        feature_sets = per_virus_data['feature_set'].unique()
        feature_mapping = create_unified_feature_mapping(feature_sets)
        
        per_virus_data['fs_short'] = per_virus_data['feature_set'].map(feature_mapping)
        if not pooled_data.empty:
            pooled_data['fs_short'] = pooled_data['feature_set'].map(feature_mapping)
            pooled_data['weight_label'] = pooled_data['class_weight'].apply(parse_weight)
        
        candidate_metrics = {
            'accuracy': 'Overall Accuracy',
            'auc': 'ROC-AUC Score',
            'recall_human': 'Sensitivity (Human-infecting)',
            'f1_human': 'F1-Score (Human-infecting)'
        }
        available_metrics = {k: v for k, v in candidate_metrics.items() if _metric_available(per_virus_data, k)}
        if not available_metrics:
            print("No valid numeric metrics found"); return False

        n_metrics = len(available_metrics)
        has_pooled = not pooled_data.empty
        
        # Widen canvas; adjust layout to improve legend placement.
        fig_width = max(16, 5 * n_metrics)
        fig, axes = plt.subplots(2, n_metrics, figsize=(fig_width, 10))
        if n_metrics == 1:
            axes = axes.reshape(-1, 1)
        
        # Define a unified color map.
        weight_colors = {
            'None': '#FF6B6B',       
            'Balanced': '#4ECDC4',   
            '5:1': '#45B7D1'          
        }
        
        for i, (metric_key, metric_title) in enumerate(available_metrics.items()):
            # Row 1: by-feature-set visualisation.
            ax1 = axes[0, i]
            try:
                plot_data = per_virus_data[['fs_short', metric_key, 'weight_label']].dropna()
                if len(plot_data) > 0:
                    unique_values = plot_data[metric_key].nunique()
                    
                    # **Change: use box plots consistently; drop violin plots.**
                    if unique_values > 5:
                        # Use box plots.
                        sns.boxplot(
                            data=plot_data,
                            x='fs_short',
                            y=metric_key,
                            hue='weight_label',
                            ax=ax1,
                            palette=weight_colors
                        )
                    else:
                        # Bar-chart fallback.
                        sns.barplot(
                            data=plot_data,
                            x='fs_short',
                            y=metric_key,
                            hue='weight_label',
                            ax=ax1,
                            palette=weight_colors,
                            alpha=0.7
                        )
                
                # **Change: use academic-style titles.**
                ax1.set_title(metric_title, fontsize=12)
                ax1.set_xlabel('Feature Set', fontsize=11)
                ax1.tick_params(axis='x', rotation=45, labelsize=10)
                
                if i == 0:
                    ax1.set_ylabel("Performance Score", fontsize=11)
                else:
                    ax1.set_ylabel("")  
                # Remove per-axes legends; add a unified legend later.
                if ax1.get_legend():
                    ax1.get_legend().remove()
                    
            except Exception as e:
                print(f"Feature plot error for {metric_key}: {e}")
                ax1.text(0.5, 0.5, f'Error: {metric_key}', ha='center', va='center', transform=ax1.transAxes)
            # Row 2: Per-virus vs pooled comparison.
            ax2 = axes[1, i]
            if has_pooled:
                try:
                    comparison_data = []
                    # Ensure all feature sets included.
                    all_features = sorted(set(per_virus_data['fs_short'].unique()) | 
                                        set(pooled_data['fs_short'].unique()))
                    main_weight = 'None'  # Select main class-weight settings to compare.
                    
                    for fs_short in all_features:
                        # Per-virus mean.
                        pv_subset = per_virus_data[
                            (per_virus_data['fs_short'] == fs_short) & 
                            (per_virus_data['weight_label'] == main_weight)
                        ]
                        if not pv_subset.empty:
                            comparison_data.append({
                                'feature': fs_short,
                                'value': pv_subset[metric_key].mean(),
                                'type': 'Per-virus'
                            })
                        
                        # Pooled results.
                        pooled_subset = pooled_data[
                            (pooled_data['fs_short'] == fs_short) & 
                            (pooled_data['weight_label'] == main_weight) &
                            (pooled_data['segment'] == 'all')
                        ]
                        if not pooled_subset.empty and metric_key in pooled_subset.columns:
                            value = pooled_subset.iloc[0][metric_key]
                            if pd.notna(value):
                                comparison_data.append({
                                    'feature': fs_short,
                                    'value': value,
                                    'type': 'Pooled'
                                })
                    
                    if comparison_data:
                        comp_df = pd.DataFrame(comparison_data)
                        
                        sns.barplot(
                            data=comp_df,
                            x='feature',
                            y='value',
                            hue='type',
                            ax=ax2,
                            palette=['#FF6B6B', '#4ECDC4'],  
                            alpha=0.8
                        )
                        # Add value labels.
                        for container in ax2.containers:
                            ax2.bar_label(container, fmt='%.2f', fontsize=8, padding=2)
                        ax2.set_title(f'{metric_title}\nPer-virus vs Pooled Evaluation', fontsize=12)
                        ax2.set_xlabel('Feature Set', fontsize=11)
                        ax2.tick_params(axis='x', rotation=45, labelsize=10)
                
                        if i == 0:
                            ax2.set_ylabel("Performance Score", fontsize=11)
                        else:
                            ax2.set_ylabel("")
                        # Remove per-axes legend.
                        if ax2.get_legend():
                            ax2.get_legend().remove()
                    else:
                        ax2.text(0.5, 0.5, 'No comparison data', ha='center', va='center', transform=ax2.transAxes)
                        
                except Exception as e:
                    print(f"Pooled comparison error for {metric_key}: {e}")
                    ax2.text(0.5, 0.5, 'Comparison error', ha='center', va='center', transform=ax2.transAxes)
            else:
                ax2.text(0.5, 0.5, 'No pooled data available', ha='center', va='center', transform=ax2.transAxes)
        
        plt.tight_layout(rect=[0, 0.05, 0.75, 0.95])
        plt.suptitle('LOVO Cross-Validation: Feature Performance and Evaluation Method Comparison', 
                     fontsize=16, y=0.97)
        # Class-weight legend.
        weight_legend_elements = []
        for weight, color in weight_colors.items():
            if weight in per_virus_data['weight_label'].values:
                weight_legend_elements.append(
                    plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.7, label=f'Class Weight: {weight}')
                )
        
        if weight_legend_elements:
            weight_legend = fig.legend(handles=weight_legend_elements, 
                                     title='Class Weight', 
                                     loc='center right', 
                                     bbox_to_anchor=(0.95, 0.75),
                                     fontsize=10)
        # Evaluation-method legend.
        if has_pooled:
            eval_legend_elements = [
                plt.Rectangle((0, 0), 1, 1, facecolor='#FF6B6B', alpha=0.8, label='Per-virus'),
                plt.Rectangle((0, 0), 1, 1, facecolor='#4ECDC4', alpha=0.8, label='Pooled')
            ]
            eval_legend = fig.legend(handles=eval_legend_elements,
                                   title='Evaluation Method',
                                   loc='center right',
                                   bbox_to_anchor=(0.95, 0.45),
                                   fontsize=10)
        # Feature-set legend.
        legend_lines = []
        for original_fs, short_fs in feature_mapping.items():
            if short_fs != original_fs:
                legend_lines.append(f"{short_fs}: {original_fs}")
        
        if legend_lines:
            legend_text = "Feature Sets:\n" + "\n".join(legend_lines[:6])
            if len(legend_lines) > 6:
                legend_text += f"\n... +{len(legend_lines)-6} more"
            
            fig.text(0.77, 0.20, legend_text, fontsize=9, verticalalignment='top',
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9))
        
        # Save.
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        
        if os.path.exists(output_file):
            print(f"✓ Academic LOVO boxplot saved: {output_file}")
            return True
        else:
            return False
            
    except Exception as e:
        print(f"✗ Academic LOVO boxplot failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def plot_simple_lovo_summary(lovo_df, output_file):
    """
    Simplified LOVO summary figure — show per-virus means vs pooled metrics
    """
    try:
        per_virus_data = lovo_df[lovo_df['result_type'] == 'per_virus'] if 'result_type' in lovo_df.columns else lovo_df
        pooled_data = lovo_df[lovo_df['result_type'] == 'pooled'] if 'result_type' in lovo_df.columns else pd.DataFrame()
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Left: distribution of human-class recall.
        if 'recall_human' in per_virus_data.columns:
            ax1 = axes[0]
            recall_data = per_virus_data['recall_human'].dropna()
            ax1.hist(recall_data, bins=20, alpha=0.7, edgecolor='black', color='lightblue')
            ax1.set_title('Per-virus Human Detection Rate Distribution')
            ax1.set_xlabel('Sensitivity (Human-infecting)')
            ax1.set_ylabel('Frequency')
            ax1.axvline(recall_data.mean(), color='red', linestyle='--', label=f'Mean: {recall_data.mean():.3f}')
            ax1.legend()
        
        # Right: per-virus vs pooled comparison (if pooled data exists).
        ax2 = axes[1]
        if not pooled_data.empty and 'recall_human' in pooled_data.columns:
            pv_mean = per_virus_data['recall_human'].mean()
            pooled_mean = pooled_data['recall_human'].mean()
            
            ax2.bar(['Per-Virus\nAverage', 'Pooled\nEvaluation'], [pv_mean, pooled_mean], 
                   color=['lightblue', 'orange'], alpha=0.8)
            ax2.set_title('Evaluation Method Comparison')
            ax2.set_ylabel('Sensitivity (Human-infecting)')
            
            # Add value labels.
            ax2.text(0, pv_mean + 0.01, f'{pv_mean:.3f}', ha='center', va='bottom')
            ax2.text(1, pooled_mean + 0.01, f'{pooled_mean:.3f}', ha='center', va='bottom')
            # Show difference.
            difference = pooled_mean - pv_mean
            ax2.text(0.5, max(pv_mean, pooled_mean) * 0.5, f'Difference:\n{difference:+.3f}', 
                    ha='center', va='center', fontsize=12, 
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))
        else:
            ax2.text(0.5, 0.5, 'No Pooled Data\nAvailable', ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Pooled Evaluation (Not Available)')
        
        plt.tight_layout()
        plt.suptitle('LOVO Evaluation: Per-virus vs Pooled Methods', fontsize=14, y=1.02)
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✓ Simple LOVO summary saved: {output_file}")
        return True
        
    except Exception as e:
        print(f"✗ Simple summary failed: {e}")
        return False

def plot_lovo_for_virology_research(lovo_df, output_file):
    """
    Simplified LOVO visualisation for virology — focus on class-weight and genome-segment contrasts.
    """
    print(f"\n--- Focused Virology LOVO Visualization ---")
    
    if lovo_df is None or lovo_df.empty:
        print("Empty LOVO data")
        return False
    
    if 'result_type' in lovo_df.columns:
        result_counts = lovo_df['result_type'].value_counts()
        print(f"Input data - Result types: {result_counts.to_dict()}")
    
    try:
        def parse_weight(weight):
            if pd.isna(weight) or weight in [None, 'None', 'nan', 'NaN']:
                return 'None'
            sw = str(weight).lower()
            if sw == 'balanced':
                return 'Balanced'
            if "'5'" in sw or '5' in sw:
                return '5:1'
            return str(weight)
        
        lovo_data = lovo_df.copy()
        lovo_data['weight_label'] = lovo_data['class_weight'].apply(parse_weight)
        
        # Split data.
        if 'result_type' in lovo_data.columns:
            per_virus_data = lovo_data[lovo_data['result_type'] == 'per_virus'].copy()
            pooled_data = lovo_data[lovo_data['result_type'] == 'pooled'].copy()
        else:
            per_virus_data = lovo_data.copy()
            pooled_data = pd.DataFrame()
        
        print(f"Separated - Per-virus: {len(per_virus_data)}, Pooled: {len(pooled_data)}")
        
        if per_virus_data.empty:
            print("No per-virus data available")
            return False
        
        # Use unified feature-set label mapping.
        def create_unified_feature_mapping(feature_sets):
            feature_mapping = {}
            for fs in sorted(feature_sets):
                if 'mixed_3mer_rscu_ctriad' in fs.lower():
                    feature_mapping[fs] = 'Mix3'
                elif 'mixed_3mer_rscu' in fs.lower():
                    feature_mapping[fs] = 'Mix2'  
                elif 'all_protein' in fs.lower():
                    feature_mapping[fs] = 'AllProt'
                elif 'ctriad' in fs.lower():
                    feature_mapping[fs] = 'CTriad'
                elif '3mer' == fs.lower():
                    feature_mapping[fs] = '3mer'
                elif '6mer' == fs.lower():
                    feature_mapping[fs] = '6mer'
                elif 'rscu' == fs.lower():
                    feature_mapping[fs] = 'RSCU'
                else:
                    feature_mapping[fs] = fs[:6] if len(fs) > 6 else fs
            return feature_mapping
        
        feature_sets = per_virus_data['feature_set'].unique()
        feature_mapping = create_unified_feature_mapping(feature_sets)
        
        per_virus_data['fs_short'] = per_virus_data['feature_set'].map(feature_mapping)
        if not pooled_data.empty:
            pooled_data['fs_short'] = pooled_data['feature_set'].map(feature_mapping)
            pooled_data['weight_label'] = pooled_data['class_weight'].apply(parse_weight)
        
        metric_definitions = {
            'accuracy':'Overall Accuracy', 'recall_human':'Sensitivity (Human-infecting)',
            'f1_human':'F1-Score (Human-infecting)', 'auc':'ROC-AUC Score'
        }
        available_metrics = [m for m in metric_definitions if _metric_available(per_virus_data, m)]
        main_metrics = [m for m in ['recall_human','auc'] if m in available_metrics][:2]
        if not main_metrics:
            print("No suitable metrics found"); return False

        main_metrics = ['recall_human', 'auc']
        main_available = [m for m in main_metrics if m in available_metrics][:2]  # At most two.
        
        fig_width = 12
        fig, axes = plt.subplots(2, 2, figsize=(fig_width, 10))
        # Define colors.
        weight_colors = {
            'None': '#FF6B6B',        
            'Balanced': '#4ECDC4',   
            '5:1': '#45B7D1'          
        }  
        # Row 1: class-weight comparison — focus.
        for i, metric in enumerate(main_available):
            metric_title = metric_definitions.get(metric, metric.replace('_', ' ').title())
            
            ax1 = axes[0, i]
            try:
                sns.boxplot(
                    data=per_virus_data,
                    x='fs_short',
                    y=metric,
                    hue='weight_label',
                    ax=ax1,
                    palette=weight_colors
                )
                
                ax1.set_title(f'{metric_title}\nClass Weight Comparison', fontsize=12)
                ax1.set_xlabel('Feature Set', fontsize=11)
                ax1.tick_params(axis='x', rotation=30, labelsize=10)
                
                if i == 0:
                    ax1.set_ylabel("Performance Score", fontsize=11)
                else:
                    ax1.set_ylabel("")
                
                # Remove legends; handle globally later.
                if ax1.get_legend():
                    ax1.get_legend().remove()
                    
            except Exception as e:
                print(f"Error in class weight plot for {metric}: {e}")
                ax1.text(0.5, 0.5, f'Error: {metric}', ha='center', va='center', transform=ax1.transAxes)
        
        # Row 2: genome-segment comparison — focus.
        # Check segment availability.
        if 'segment' in per_virus_data.columns:
            available_segments = per_virus_data['segment'].unique()
            print(f"Available segments for analysis: {available_segments}")
            
            for i, metric in enumerate(main_available):
                ax2 = axes[1, i]
                try:
                    # Use segments with sufficient data only.
                    segment_data = []
                    for segment in available_segments:
                        segment_subset = per_virus_data[per_virus_data['segment'] == segment]
                        if len(segment_subset) >= 5:  # Require at least 5 data points.
                            segment_data.append(segment)
                    
                    if len(segment_data) >= 2:  # Need at least two segments to compare.
                        plot_data = per_virus_data[per_virus_data['segment'].isin(segment_data)]
                        sns.boxplot(
                            data=plot_data,
                            x='segment',
                            y=metric,
                            hue='weight_label',
                            ax=ax2,
                            palette=weight_colors
                        )
                        
                        metric_title = metric_definitions.get(metric, metric.replace('_', ' ').title())
                        ax2.set_title(f'{metric_title}\nGenome Segment Comparison', fontsize=12)
                        ax2.set_xlabel('Genome Segment', fontsize=11)
                        ax2.tick_params(axis='x', rotation=0, labelsize=10)
                        
                        if i == 0:
                            ax2.set_ylabel("Performance Score", fontsize=11)
                        else:
                            ax2.set_ylabel("")
                        # Remove legend.
                        if ax2.get_legend():
                            ax2.get_legend().remove()
                        
                        print(f"✓ Segment analysis for {metric}: {len(segment_data)} segments")
                    else:
                        ax2.text(0.5, 0.5, f'Insufficient segment data\n({len(segment_data)} segments with enough data)', 
                                ha='center', va='center', transform=ax2.transAxes)
                        ax2.set_title('Segment Analysis\n(Insufficient Data)', fontsize=12)
                        
                except Exception as e:
                    print(f"Error in segment analysis for {metric}: {e}")
                    ax2.text(0.5, 0.5, f'Segment analysis error:\n{str(e)[:50]}...', 
                            ha='center', va='center', transform=ax2.transAxes)
        else:
            # If no segment info.
            for i in range(2):
                axes[1, i].text(0.5, 0.5, 'No segment information\navailable in data', 
                               ha='center', va='center', transform=axes[1, i].transAxes)
                axes[1, i].set_title('Segment Analysis\n(No Data)', fontsize=12)
        
        # Adjust layout.
        plt.tight_layout(rect=[0, 0.03, 0.85, 0.96])
        plt.suptitle('LOVO Analysis: Class Weight and Genome Segment Effects', 
                     fontsize=16, y=0.98)
        # Right-side legend — show class weight only.
        weight_legend_elements = []
        for weight, color in weight_colors.items():
            if weight in per_virus_data['weight_label'].values:
                weight_legend_elements.append(
                    plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.7, label=weight)
                )
        
        if weight_legend_elements:
            fig.legend(handles=weight_legend_elements,
                      title='Class Weight',
                      loc='center right',
                      bbox_to_anchor=(0.98, 0.5),
                      fontsize=10,
                      title_fontsize=11)
        # Feature-set notes.
        if feature_mapping:
            legend_lines = []
            for original_fs, short_fs in feature_mapping.items():
                if short_fs != original_fs:
                    legend_lines.append(f"{short_fs}: {original_fs}")
            if legend_lines:
                # Limit number displayed.
                display_lines = legend_lines[:4]
                legend_text = "Feature Sets:\n" + "\n".join(display_lines)
                if len(legend_lines) > 4:
                    legend_text += f"\n+{len(legend_lines)-4} more..."
                
                fig.text(0.87, 0.15, legend_text, 
                        fontsize=8, verticalalignment='top',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))
        
        # Save.
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        
        if os.path.exists(output_file):
            print(f"✓ Focused virology LOVO visualization saved: {output_file}")
            print(f" Focused Analysis Summary:")
            print(f"  - Class weight comparison: {len(main_available)} metrics")
            print(f"  - Segment comparison: {len(main_available)} metrics")
            
            return True
        else:
            return False
            
    except Exception as e:
        print(f"✗ Focused virology visualization failed: {e}")
        import traceback
        traceback.print_exc()
        return False
   
def generate_lovo_summary_table(lovo_df, output_file):
    """
    Generate a concise comparison table for LOVO results.
    """
    try:
        per_virus_data = lovo_df[lovo_df['result_type'] == 'per_virus']
        pooled_data = lovo_df[lovo_df['result_type'] == 'pooled']
        
        summary_rows = []
        
        for fs in per_virus_data['feature_set'].unique():
            for cw in per_virus_data['class_weight'].unique():
                # Per-virus mean.
                pv_subset = per_virus_data[
                    (per_virus_data['feature_set'] == fs) & 
                    (per_virus_data['class_weight'] == cw)
                ]
                # Pooled results (all segments).
                pooled_subset = pooled_data[
                    (pooled_data['feature_set'] == fs) & 
                    (pooled_data['class_weight'] == cw) &
                    (pooled_data['segment'] == 'all')
                ]
                if not pv_subset.empty:
                    row = {
                        'feature_set': fs,
                        'class_weight': cw,
                        'n_viruses': len(pv_subset),
                        'total_samples': pv_subset['n_samples'].sum(),
                    }            
                    # Per-virus metrics.
                    for metric in ['accuracy', 'auc', 'recall_human', 'recall_nonhuman']:
                        if metric in pv_subset.columns:
                            row[f'{metric}_pv_mean'] = pv_subset[metric].mean()
                            row[f'{metric}_pv_std'] = pv_subset[metric].std()
                    # Pooled metrics.
                    if not pooled_subset.empty:
                        pooled_row = pooled_subset.iloc[0]
                        for metric in ['accuracy', 'auc', 'recall_human', 'recall_nonhuman']:
                            if metric in pooled_subset.columns:
                                row[f'{metric}_pooled'] = pooled_row[metric]
                                pv_mean = row.get(f'{metric}_pv_mean', np.nan)
                                if not pd.isna(pv_mean) and not pd.isna(pooled_row[metric]):
                                    row[f'{metric}_diff'] = pooled_row[metric] - pv_mean
                    
                    summary_rows.append(row)
        
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"✓ LOVO summary table saved: {output_file}")
        return True
        
    except Exception as e:
        print(f"✗ Failed to generate summary table: {e}")
        return False

def print_lovo_key_findings(lovo_df):
    print(f"\n--- Key Findings Summary ---")
    try:
        if 'result_type' in lovo_df.columns:
            per_virus_data = lovo_df[lovo_df['result_type'] == 'per_virus']
            pooled_data = lovo_df[lovo_df['result_type'] == 'pooled']
        else:
            per_virus_data = lovo_df
            pooled_data = pd.DataFrame()
        if not per_virus_data.empty:
            # Ability to detect human-infecting viruses.
            if 'recall_human' in per_virus_data.columns:
                human_recall = per_virus_data['recall_human'].dropna()
                if len(human_recall) > 0:
                    avg_sensitivity = human_recall.mean()
                    print(f" Human Pathogen Detection Sensitivity: {avg_sensitivity:.1%}")
                    
                    # Classification performance assessment.
                    excellent = (human_recall >= 0.8).sum()
                    good = ((human_recall >= 0.6) & (human_recall < 0.8)).sum()
                    poor = (human_recall < 0.6).sum()
                    total = len(human_recall)
                    
                    print(f"   • Excellent performance (≥80%): {excellent}/{total} viruses ({excellent/total:.1%})")
                    print(f"   • Good performance (60-80%): {good}/{total} viruses ({good/total:.1%})")
                    print(f"   • Needs improvement (<60%): {poor}/{total} viruses ({poor/total:.1%})")
            
            # Overall discriminative ability.
            if 'auc' in per_virus_data.columns:
                auc_values = per_virus_data['auc'].dropna()
                if len(auc_values) > 0:
                    avg_auc = auc_values.mean()
                    print(f" Average Discrimination Ability (AUC): {avg_auc:.3f}")
                    
                    strong_discrimination = (auc_values >= 0.8).sum()
                    moderate_discrimination = ((auc_values >= 0.7) & (auc_values < 0.8)).sum()
                    
                    print(f"   • Strong discrimination (≥0.8): {strong_discrimination}/{len(auc_values)} viruses")
                    print(f"   • Moderate discrimination (0.7-0.8): {moderate_discrimination}/{len(auc_values)} viruses")
            
            # Identify the best feature sets.
            if 'feature_set' in per_virus_data.columns and 'recall_human' in per_virus_data.columns:
                fs_performance = per_virus_data.groupby('feature_set')['recall_human'].mean().sort_values(ascending=False)
                if len(fs_performance) > 0:
                    best_fs = fs_performance.index[0]
                    best_score = fs_performance.iloc[0]
                    print(f" Best Feature Set: '{best_fs}' (Human detection: {best_score:.1%})")
            
            # Pooled vs per-virus comparison.
            if not pooled_data.empty and 'recall_human' in pooled_data.columns:
                pooled_all = pooled_data[pooled_data['segment'] == 'all']
                if not pooled_all.empty:
                    pooled_sensitivity = pooled_all['recall_human'].mean()
                    per_virus_sensitivity = per_virus_data['recall_human'].mean()  
                    if not pd.isna(pooled_sensitivity) and not pd.isna(per_virus_sensitivity):
                        diff = pooled_sensitivity - per_virus_sensitivity
                        print(f" Real-world vs Lab Performance:")
                        print(f"   • Pooled (real-world): {pooled_sensitivity:.1%}")
                        print(f"   • Per-virus average: {per_virus_sensitivity:.1%}")
                        print(f"   • Difference: {diff:+.1%} ({'Better' if diff > 0 else 'Worse'} in real-world)")
            
        print(f" Main results visualized in 'lovo_virology_analysis.png'!")
        
    except Exception as e:
        print(f" Error generating key findings: {e}")

def plot_feature_importance_comparison(df, feature_sets_dict, output_prefix):
    print(f"\n--- Starting Feature Importance Analysis ---")
    
    if df is None or df.empty:
        print(" Empty dataframe for feature importance")
        return False
    
    print(f"Input data shape: {df.shape}")
    print(f"Feature sets to analyze: {list(feature_sets_dict.keys())}")
    
    try:
        from sklearn.ensemble import RandomForestClassifier
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import os
        from datetime import datetime
        # Check label columns.
        if 'final_label' not in df.columns:
            print("✗ Missing 'final_label' column")
            return False
        unique_labels = df['final_label'].unique()
        print(f"Unique labels: {unique_labels}")
        # Focus on mixed feature sets and a few main sets.
        priority_sets = ['mixed_3mer_rscu_ctriad', 'mixed_3mer_rscu', '3mer', 'rscu', 'ctriad']
        selected_sets = {}
        # Prefer feature sets in priority_sets.
        for feat_name in priority_sets:
            if feat_name in feature_sets_dict:
                selected_sets[feat_name] = feature_sets_dict[feat_name]
        # If fewer than six, add others.
        remaining_sets = {k: v for k, v in feature_sets_dict.items() if k not in selected_sets}
        for feat_name, feat_func in list(remaining_sets.items())[:6-len(selected_sets)]:
            selected_sets[feat_name] = feat_func
        print(f"Selected feature sets for analysis: {list(selected_sets.keys())}")
        
        top_features_summary = []
        n_sets = len(selected_sets)
        
        # Create table-like layout, one row per feature set.
        fig, axes = plt.subplots(n_sets, 1, figsize=(16, 3.5 * n_sets))
        if n_sets == 1:
            axes = [axes]
        
        # Define colors by feature type (for mixed sets).
        feature_type_colors = {
            '3mer': '#FF6B6B',      
            'rscu': '#4ECDC4',       
            'protein': '#45B7D1',   
            'ctriad': '#FFA500',    
            'other': '#999999'      
        }
        
        def classify_feature_type(feature_name):
            """Infer feature type from feature names."""
            feature_lower = feature_name.lower()
            if any(x in feature_lower for x in ['3mer', 'trigram', 'codon']):
                return '3mer'
            elif 'rscu' in feature_lower:
                return 'rscu'
            elif any(x in feature_lower for x in ['protein', 'aa', 'amino', 'ctriad']):
                if 'ctriad' in feature_lower:
                    return 'ctriad'
                else:
                    return 'protein'
            else:
                return 'other'
        
        success_count = 0
        has_mixed_sets = False  # Track whether mixed sets exist.
        used_feature_types = set()  # Track which feature types are used.
        
        for idx, (feat_name, feat_func) in enumerate(selected_sets.items()):
            print(f"\nProcessing feature set: {feat_name}")
            
            try:
                # Get feature columns.
                if callable(feat_func):
                    feat_cols = feat_func(df)
                else:
                    feat_cols = feat_func  # If already a list of column names.
                
                print(f"Feature columns count: {len(feat_cols)}")
                
                if len(feat_cols) < 5:
                    print(f"Skipping {feat_name}: insufficient features ({len(feat_cols)})")
                    axes[idx].text(0.5, 0.5, f'{feat_name}\n(insufficient features)', 
                                  ha='center', va='center', transform=axes[idx].transAxes)
                    axes[idx].set_xticks([])
                    axes[idx].set_yticks([])
                    continue
                
                # Check feature columns exist.
                missing_cols = [col for col in feat_cols if col not in df.columns]
                if missing_cols:
                    print(f"Warning: {len(missing_cols)} missing columns in {feat_name}")
                    feat_cols = [col for col in feat_cols if col in df.columns]
                
                if len(feat_cols) < 5:
                    print(f"Skipping {feat_name}: too few valid features after filtering")
                    continue
                
                # Prepare data.
                X = df[feat_cols].fillna(0).values
                y = (df['final_label'] == 'human').astype(int)
                
                print(f"Training data shape: X={X.shape}, y={y.shape}")
                print(f"Label distribution: {np.bincount(y)}")
                
                # Train Random Forest.
                clf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
                clf.fit(X, y)
                
                # Get feature importances.
                importances = clf.feature_importances_
                indices = np.argsort(importances)[::-1][:10]  # Take top 10 only.
                
                print(f"Top importance score: {importances[indices[0]]:.4f}")
                
                # Assign colors for mixed sets and record used types.
                is_mixed = 'mixed' in feat_name.lower()
                if is_mixed:
                    has_mixed_sets = True
                    colors = []
                    current_set_types = set()
                    for i in indices:
                        feature_type = classify_feature_type(feat_cols[i])
                        colors.append(feature_type_colors[feature_type])
                        used_feature_types.add(feature_type)
                        current_set_types.add(feature_type)
                else:
                    # Use a single color for non-mixed sets.
                    colors = ['#666666'] * len(indices)
                    current_set_types = set()
                
                # Save top-feature info.
                for rank, i in enumerate(indices):
                    feature_type = classify_feature_type(feat_cols[i]) if is_mixed else 'single_type'
                    top_features_summary.append({
                        'feature_set': feat_name,
                        'rank': rank + 1,
                        'feature': feat_cols[i], 
                        'importance': importances[i],
                        'feature_type': feature_type
                    })
                
                # Plot horizontal bar charts to show feature names.
                ax = axes[idx]
                y_pos = np.arange(len(indices))
                bars = ax.barh(y_pos, importances[indices], color=colors, alpha=0.8)
                # Set y-tick labels to feature names (simplified).
                feature_labels = []
                for i in indices:
                    label = feat_cols[i]
                    if len(label) > 30:  # Truncate overly long names.
                        label = label[:27] + "..."
                    feature_labels.append(label)
                
                ax.set_yticks(y_pos)
                ax.set_yticklabels(feature_labels, fontsize=9)
                ax.set_xlabel('Feature Importance', fontsize=11)
                
                short_name = feat_name.replace('mixed_', '').replace('_', '-')
                ax.set_title(f'Top 10 Features: {short_name}', fontsize=12, fontweight='bold')
                
                # Add type counts inside subplot for mixed sets.
                if is_mixed and current_set_types:
                    # Count features by type.
                    type_counts = {}
                    for i in indices:
                        ftype = classify_feature_type(feat_cols[i])
                        type_counts[ftype] = type_counts.get(ftype, 0) + 1
                    
                    # Show type counts in the subplot’s top-right.
                    legend_text = "Types:\n"
                    for ftype in sorted(current_set_types):
                        count = type_counts.get(ftype, 0)
                        if count > 0:
                            color_hex = feature_type_colors[ftype]
                            legend_text += f"● {ftype.upper()}: {count}\n"
                    
                    ax.text(0.98, 0.98, legend_text.strip(), 
                           transform=ax.transAxes, fontsize=9,
                           verticalalignment='top', horizontalalignment='right',
                           bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.9, edgecolor='gray'))
                
                success_count += 1
                print(f"✓ Feature importance plot created for {feat_name}")
                
            except Exception as e:
                print(f"✗ Error processing feature set {feat_name}: {e}")
                axes[idx].text(0.5, 0.5, f'{feat_name}\n(error)', 
                              ha='center', va='center', transform=axes[idx].transAxes)
                axes[idx].set_xticks([])
                axes[idx].set_yticks([])
        
        # Tighten layout; reserve space for right-side legend.
        plt.tight_layout(rect=[0, 0.05, 0.82, 0.95])
        
        # Add a global title. 
        fig.suptitle('Feature Importance Analysis: Top 10 Features by Model Type', 
                     fontsize=16, fontweight='bold', y=0.98)
        
        # Add global legend only when mixed sets exist.
        if has_mixed_sets and used_feature_types:
            # Create global legend for the types actually used.
            legend_elements = []
            for ftype in sorted(used_feature_types):
                if ftype != 'other':  #  Usually omit the “other” type.
                    color = feature_type_colors[ftype]
                    label = ftype.upper()
                    legend_elements.append(
                        plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.8, label=label)
                    )
            if legend_elements:  # Create legend only if there are entries.
                legend = fig.legend(handles=legend_elements, 
                          title='Feature Types\n(Mixed Models Only)',
                          loc='center right', 
                          bbox_to_anchor=(0.98, 0.5),
                          fontsize=10,
                          title_fontsize=11,
                          frameon=True,
                          fancybox=True,
                          shadow=True)
                legend.get_frame().set_facecolor('white')
                legend.get_frame().set_alpha(0.9)    
        # Add methods note.
        method_text = ("Method: Random Forest (n=100)\n"
                      "Ranking: Feature Importance Score\n"
                      "Colors indicate feature types in mixed models")
        
        fig.text(0.84, 0.15, method_text, fontsize=9, 
                verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8))
        # Save figure.
        plot_file = f'{output_prefix}_top10_features.png'
        os.makedirs(os.path.dirname(plot_file), exist_ok=True)
        plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        # Verify image file creation.
        plot_success = False
        if os.path.exists(plot_file):
            file_size = os.path.getsize(plot_file)
            timestamp = datetime.fromtimestamp(os.path.getmtime(plot_file))
            print(f"✓ Top 10 feature importance plot saved: {plot_file}")
            print(f"  File size: {file_size} bytes")
            print(f"  Timestamp: {timestamp}")
            print(f"  Mixed sets processed: {has_mixed_sets}")
            print(f"  Feature types found: {sorted(used_feature_types) if used_feature_types else 'None'}")
            plot_success = True
        else:
            print(f"✗ Feature importance plot file not created: {plot_file}")
        
        # Save as a supplementary table.
        csv_success = False
        if top_features_summary:
            # Create a detailed supplementary table.
            summary_df = pd.DataFrame(top_features_summary)
            # Sort by feature set and rank.
            summary_df = summary_df.sort_values(['feature_set', 'rank'])
            # Format importance scores.
            summary_df['importance_formatted'] = summary_df['importance'].apply(lambda x: f"{x:.4f}")
            # Reorder columns.
            columns_order = ['feature_set', 'rank', 'feature', 'importance', 'importance_formatted', 'feature_type']
            summary_df = summary_df[columns_order]
            
            csv_file = f"{output_prefix}_supplementary_table.csv"
            summary_df.to_csv(csv_file, index=False)
            
            if os.path.exists(csv_file):
                file_size = os.path.getsize(csv_file)
                timestamp = datetime.fromtimestamp(os.path.getmtime(csv_file))
                print(f"✓ Supplementary feature table saved: {csv_file}")
                print(f"  File size: {file_size} bytes")
                print(f"  Timestamp: {timestamp}")
                print(f"  Total features recorded: {len(summary_df)}")
                csv_success = True
            else:
                print(f"✗ Supplementary table CSV not created: {csv_file}")
        
        return plot_success or csv_success
        
    except Exception as e:
        print(f"✗ Feature importance analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False
       
def plot_clustermap(
    data, title, outpath, figsize=(22, 12), cmap='magma',
    font_scale=0.7, dendrogram_ratio=(.15, .08),
    col_cluster=True, col_linkage=None
):
    # Data cleaning.
    data = data.loc[~data.isnull().all(axis=1), :]    
    data = data.loc[:, ~data.isnull().all(axis=0)]    
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(axis=0, how='any')             
    data = data.dropna(axis=1, how='any')             
    
    if data.shape[1] == 0 or data.shape[0] == 0:
        print(f"[Warning] Empty data for {title}, skip.")
        return None
    if data is None or data.empty:
        print(f"[Warning] Empty data for {title}, skip.")
        return None
    # Guard against all-NaN/all-Inf.
    data = data.replace([np.inf, -np.inf], np.nan)
    if data.isnull().all().all():
        print(f"[Warning] All values NaN/Inf in {title}, skip.")
        return None
    # Standardise data.
    row_means = data.mean(axis=1)
    row_stds = data.std(axis=1).replace(0, 1)
    data_norm = data.sub(row_means, axis=0).div(row_stds, axis=0)
    sns.set_theme(font_scale=font_scale)
    
    # Create figure.
    plt.figure(figsize=figsize)
    plt.clf()
    # Create clustermap with automatic colorbar disabled.
    g = sns.clustermap(
        data_norm,
        figsize=figsize,
        cmap=cmap,
        cbar=False,  # Disable the auto colorbar.
        col_cluster=col_cluster,
        col_linkage=col_linkage,
        row_cluster=True,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=dendrogram_ratio,
        colors_ratio=0.03,
    )
    
    # Add colorbar manually.
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize
    
    heatmap_bbox = g.ax_heatmap.get_position()
    new_left = 0.20
    new_width = 0.55
    new_height = heatmap_bbox.height
    new_bottom = heatmap_bbox.y0
    
    g.ax_heatmap.set_position([new_left, new_bottom, new_width, new_height])
    
    # Adjust dendrogram position.
    if g.ax_row_dendrogram is not None:
        dendro_bbox = g.ax_row_dendrogram.get_position()
        dendro_width = dendro_bbox.width
        dendro_left = new_left - dendro_width - 0.005
        g.ax_row_dendrogram.set_position([dendro_left, dendro_bbox.y0, dendro_width, dendro_bbox.height])
    
    if g.ax_col_dendrogram is not None:
        col_dendro_bbox = g.ax_col_dendrogram.get_position()
        g.ax_col_dendrogram.set_position([new_left, col_dendro_bbox.y0, new_width, col_dendro_bbox.height])
    
    # Get final heatmap position and create colorbar.
    final_heatmap_pos = g.ax_heatmap.get_position()
    cbar_height = final_heatmap_pos.height * 0.8
    cbar_y = final_heatmap_pos.y0 + (final_heatmap_pos.height - cbar_height) / 2
    cbar_left = final_heatmap_pos.x0 + final_heatmap_pos.width + 0.08
    cbar_ax = g.figure.add_axes([cbar_left, cbar_y, 0.02, cbar_height])
    
    # Create colorbar.
    vmin, vmax = data_norm.min().min(), data_norm.max().max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    cbar = ColorbarBase(cbar_ax, cmap=plt.cm.get_cmap(cmap), norm=norm, orientation='vertical')
    cbar.set_label('Normalized Expression', rotation=270, labelpad=15, fontsize=9)
    
    # Set labels and titles.
    g.ax_heatmap.set_xlabel("Top Variable Features", fontsize=12)
    g.ax_heatmap.set_ylabel("Virus Species", fontsize=12)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=7)
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), 
        rotation=0, fontsize=6, ha='left', va='center'
    )
    g.ax_heatmap.tick_params(axis='y', which='major', pad=5)
    
    # Set title.
    g.figure.suptitle(title, y=0.95, fontsize=16)
    
    # Hide extra axes.
    for ax in g.figure.get_axes():
        if ax not in [g.ax_heatmap, g.ax_row_dendrogram, g.ax_col_dendrogram, cbar_ax]:
            ax.set_visible(False)
    
    # Save image.
    g.savefig(outpath, dpi=220, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    print(f"[Info] Saved: {outpath}")
    return g

def top_var_features(df, topn=30):
    """Select features with highest variance."""
    if df.shape[1] <= topn:
        return df
    var_cols = df.var().sort_values(ascending=False).head(topn).index
    return df[var_cols]

def main_visualization_section(merged_features, train_test_df, lovo_df, feature_sets): 
    print("\n==================================================")
    print("GENERATING VISUALIZATIONS")
    print("==================================================")
    
    # 0) Environment checks.
    print("\n--- Environment Check ---")
    if not debug_plot_environment():
        print("✗ Plotting environment check failed! Cannot proceed with visualizations.")
        return False
    
    # Prepare dataset for visualisation.
    df_for_viz = merged_features['intersect']
    print(f"Using intersect dataset for visualizations ({len(df_for_viz)} samples)")
    
    # Data validation.
    print(f"Data shape: {df_for_viz.shape}")
    print(f"Columns sample: {df_for_viz.columns.tolist()[:10]}...")  # Show the first 10 columns.
    if 'final_label' in df_for_viz.columns:
        print(f"Labels: {df_for_viz['final_label'].value_counts().to_dict()}")
    else:
        print(" Warning: 'final_label' column not found!")
    
    # Create results directory.
    os.makedirs('results', exist_ok=True)
    
    success_count = 0
    total_operations = 0
    
    # --- 1) Performance comparison line plot ---
    print("\n--- Performance Comparison Plots ---")
    if not train_test_df.empty:
        total_operations += 1
        if plot_performance_comparison_single(train_test_df, "results/performance_comparison_single"):
            success_count += 1
            print("✓ Performance comparison plots generated")
        else:
            print("✗ Performance comparison plots failed")
    else:
        print(" No train-test results available for performance plots")
    
    # --- 2) Visualise LOVO results ---
    print("\n--- LOVO Performance Plots ---")
    if lovo_df is not None and not lovo_df.empty:
        total_operations += 1
        if plot_lovo_boxplot(lovo_df, "results/lovo_performance_boxplot.png"):
            success_count += 1
            print("✓ LOVO boxplot generated")
        else:
            print("✗ LOVO boxplot failed")
    else:
        print(" No LOVO results available for visualization")
    
    # --- 3) UMAP projection ---
    print("\n--- UMAP Dimensionality Reduction ---")
    
    umap_configs = [
        ('3mer_', "UMAP of 3-mer Features", "results/umap_3mer_features.png"),
        ('6mer_', "UMAP of 6-mer Features", "results/umap_6mer_features.png"),
        ('RSCU_', "UMAP of RSCU Features", "results/umap_rscu_features.png"),
        (('CTDC', 'CTDT', 'CTDD', 'CTriad', 'DistancePair', 'PseAAC'), 
         "UMAP of Protein Features", "results/umap_protein_features.png"),
        (('3mer_', 'RSCU_', 'CTriad'), 
         "UMAP of Mixed Features (3-mer + RSCU + CTriad)", "results/umap_mixed_features.png")
    ]
    
    for prefix, title, output_file in umap_configs:
        total_operations += 1
        if isinstance(prefix, tuple):
            # Multiple-prefix case.
            cols = [c for c in df_for_viz.columns if any(c.startswith(p) for p in prefix)]
        else:
            # Single-prefix case.
            cols = [c for c in df_for_viz.columns if c.startswith(prefix)]
        
        if len(cols) > 0:
            if plot_umap(df_for_viz, cols, title, output_file):
                success_count += 1
        else:
            print(f" No columns found for {title}")
    
    print("✓ UMAP plots completed")
    
    # --- 4) Feature-importance analysis ---
    print("\n--- Feature Importance Analysis ---")
    total_operations += 1
    if plot_feature_importance_comparison(df_for_viz, feature_sets, "results/feature_importance"):
        success_count += 1
        print("✓ Feature importance plots completed")
    else:
        print("✗ Feature importance plots failed")
    
    # Summary.
    print(f"\n--- Visualization Summary ---")
    print(f"Total operations: {total_operations}")
    print(f"Successful operations: {success_count}")
    print(f"Success rate: {success_count/total_operations*100:.1f}%" if total_operations > 0 else "N/A")
    
    if success_count == 0:
        print("✗ All visualization operations failed!")
        return False
    elif success_count < total_operations:
        print(f" Some visualization operations failed ({total_operations - success_count} failures)")
        return True
    else:
        print("✓ All visualization operations completed successfully!")
        return True

# --- Main pipeline function ---
def run_pipeline(config):
    """Run the full analysis pipeline."""
    # Concatenate human/nonhuman FASTA into a single file.
    print("Combining FASTA files for feature extraction...")
    combined_fasta_path = 'sequences_combined.fasta'
    with open(combined_fasta_path, 'w') as outfile:
        with open(config['human_fasta']) as infile:
            outfile.write(infile.read())
        with open(config['nonhuman_fasta']) as infile:
            outfile.write(infile.read())
    print(f"✓ Combined FASTA created: {combined_fasta_path}")
    
    # ========== Phase 1: build main label table ==========
    meta_main = build_meta_main(
        config['human_fasta'],
        config['nonhuman_fasta'],
        'meta_main_post_relabel.csv'
    )
    
    # ========== Phase 2: FASTA preprocessing & alignment ==========
    # Process CDS.
    if 'cds_fasta' in config:
        cds_meta = label_cds_fasta(
            config['cds_fasta'],
            meta_main,
            'cds_meta_labeled.csv'
        )
    
    # Process protein.
    if 'protein_fasta' in config:
        protein_meta = label_protein_fasta(
            config['protein_fasta'],
            meta_main,
            'protein_meta_labeled.csv'
        )
    
    # ========== Phase 3: feature engineering ==========
    print("\n=== Phase 3: Feature Engineering and Machine Learning ===")
    
    # Extract k-mer features.
    kmer_df = extract_kmer_features(
        combined_fasta_path,  
        meta_main,
        k_values=[2, 3, 4, 5, 6],
        output_csv='kmer_features.csv'
    )
    
    # Extract RSCU/CAI features.
    rscu_df = None
    if 'cds_fasta' in config and config['cds_fasta']:
        rscu_df = extract_rscu_cai_features(
            cds_fasta=config['cds_fasta'],
            meta_main=meta_main,
            reference_fasta=config.get('reference_fasta', 'human_HK_CDS.cleaned.fasta'),
            output_csv='rscu_cai_features.csv',
            min_partial_length=600  # Can set to 300 to test gain.
        )
    
    # Extract protein features.
    protein_df = None
    if 'protein_fasta' in config and config['protein_fasta']:
        protein_df = extract_protein_features(
            config['protein_fasta'],
            protein_meta,
            meta_main, 
            output_csv='protein_features.csv'
        )

    meta_main['id_clean'] = meta_main['id'].str.replace(r'\.\d+$', '', regex=True)

    if rscu_df is not None and not rscu_df.empty:
        rscu_df['nucleotide_id_clean'] = rscu_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
    else:
        rscu_df = pd.DataFrame(columns=['nucleotide_id', 'nucleotide_id_clean', 'species_std', 'segment', 'final_label'])

    if protein_df is not None and not protein_df.empty and 'nucleotide_id' in protein_df.columns:
        protein_df['nucleotide_id_clean'] = protein_df['nucleotide_id'].str.replace(r'\.\d+$', '', regex=True)
    else:
        protein_df = pd.DataFrame(columns=['protein_id', 'nucleotide_id', 'nucleotide_id_clean', 'species_std', 'segment', 'final_label'])

    # Only perform alignment/merge (if required) when the corresponding tables are non-empty
    if not rscu_df.empty:
        common_rscu_ids = set(meta_main['id_clean']) & set(rscu_df['nucleotide_id_clean'])
        rscu_only = rscu_df[rscu_df['nucleotide_id_clean'].isin(common_rscu_ids)].copy()
        rscu_only = rscu_only.merge(
            meta_main[['id_clean', 'species_std']],
            left_on='nucleotide_id_clean', right_on='id_clean', how='left'
        )

    if not protein_df.empty:
        common_protein_ids = set(meta_main['id_clean']) & set(protein_df['nucleotide_id_clean'])
        protein_only = protein_df[protein_df['nucleotide_id_clean'].isin(common_protein_ids)].copy()
        protein_only = protein_only.merge(
            meta_main[['id_clean', 'species_std']],
            left_on='nucleotide_id_clean', right_on='id_clean', how='left'
        )
    # Merge all features.
    merged_features = merge_features(
        meta_main,
        kmer_df,
        (None if rscu_df.empty else rscu_df),
        (None if protein_df.empty else protein_df)
    )
    # Add pipeline_feature_set_coverage_summary.csv
    summary = {
        'Feature': [],
        'Sample_Count': [],
        'Feature_Count': [],
        'Human_Count': [],
        'Nonhuman_Count': [],
    }
    for name, df in merged_features.items():
        if df is None or df.empty: 
            continue
        summary['Feature'].append(name)
        summary['Sample_Count'].append(len(df))
        summary['Feature_Count'].append(len([c for c in df.columns if c not in ['id', 'nucleotide_id', 'species_std', 'segment', 'final_label']]))
        summary['Human_Count'].append((df['final_label'] == 'human').sum())
        summary['Nonhuman_Count'].append((df['final_label'] == 'nonhuman').sum())
    pd.DataFrame(summary).to_csv("pipeline_feature_set_coverage_summary.csv", encoding='utf-8', index=False)
    print("✓ Summary table saved: pipeline_feature_set_coverage_summary.csv")
    
    # ========== Machine-learning evaluation ==========
    print("\n" + "="*50)
    print("MACHINE LEARNING EVALUATION")
    print("="*50)
    
    # Check the intersect dataset.
    if 'intersect' not in merged_features or merged_features['intersect'] is None or merged_features['intersect'].empty:
        print("Error: 'intersect' dataset not found or empty. Cannot proceed with ML evaluation.")
        return
    
    # ========== Test: run only 12 human + 8 nonhuman species ==========

    intersect_df = merged_features['intersect']

    # Sample species.
    human_species = intersect_df[intersect_df['final_label'] == 'human']['species_std'].drop_duplicates()
    nonhuman_species = intersect_df[intersect_df['final_label'] == 'nonhuman']['species_std'].drop_duplicates()

    # Sampling config for debugging; disable by commenting out the whole if-block.
    num_human = 12
    num_nonhuman = 8

    if len(human_species) >= num_human and len(nonhuman_species) >= num_nonhuman:
        selected_species = pd.concat([
            human_species.sample(num_human, random_state=42),
            nonhuman_species.sample(num_nonhuman, random_state=42)
        ])
        print("[Debug] Selected species:", selected_species.tolist())

        # Filter intersect.
        intersect_df = intersect_df[intersect_df['species_std'].isin(selected_species)]
        merged_features['intersect'] = intersect_df

        # Quickly check label distribution.
        print("[Debug] Filtered label counts:\n", intersect_df['final_label'].value_counts())
    else:
        print(f"[Warning] Too few species to sample {num_human} human + {num_nonhuman} nonhuman. Running full set.")

    # ========== Debug End ==========

    # Run ML evaluation.
    train_test_df, lovo_df = run_ml_evaluation(
        merged_features, 
        output_dir='results',
        run_lovo=True,
        max_species=None,  # Or set numerical limits.
        optimize_params=True  # Enable group-aware 5-fold parameter tuning.
    )
    
    if train_test_df is None:
        print("Error: Machine learning evaluation failed. Exiting.")
        return
    
    print("\n==================================================")
    print("GENERATING VISUALIZATIONS") 
    print("==================================================")
    
    # Keep these variable definitions for the visualization section.
    df_for_viz = merged_features['intersect']
    print(f"Using intersect dataset for visualizations ({len(df_for_viz)} samples)")
    # Create results directory.
    os.makedirs('results', exist_ok=True)
    # 0) Environment checks.
    print("\n--- Environment Check ---")
    if not debug_plot_environment():
        print("✗ Plotting environment check failed! Cannot proceed with visualizations.")
    # Data validation.
    print(f"Data shape: {df_for_viz.shape}")
    print(f"Columns sample: {df_for_viz.columns.tolist()[:10]}...")
    if 'final_label' in df_for_viz.columns:
        print(f"Labels: {df_for_viz['final_label'].value_counts().to_dict()}")
    else:
        print(" Warning: 'final_label' column not found!")
    
    success_count = 0
    total_operations = 0
    
    # --- 1) Performance comparison line plot ---
    print("\n--- Performance Comparison Plots ---")
    if not train_test_df.empty:
        total_operations += 1
        if plot_performance_comparison_single(train_test_df, "results/performance_comparison_single"):
            print("✓ Performance comparison plots generated")
        else:
            print("✗ Performance comparison plots failed")
    else:
        print(" No train-test results available for performance plots")
    
    # --- 2) Visualise LOVO results ---
    print("\n--- LOVO Performance Analysis ---")
    
    # Read full data directly from CSV.
    try:
        complete_lovo_df = pd.read_csv('results/lovo_results.csv')
        print(f" Loaded complete LOVO data from CSV: {complete_lovo_df.shape}")
        
        if 'result_type' in complete_lovo_df.columns:
            type_counts = complete_lovo_df['result_type'].value_counts(dropna=False)
            print(f"Result type distribution: {type_counts.to_dict()}")
        
        lovo_data_for_plotting = complete_lovo_df
        
    except Exception as e:
        print(f" Failed to read CSV: {e}")
        lovo_data_for_plotting = lovo_df
    
    if lovo_data_for_plotting is not None and not lovo_data_for_plotting.empty:
        total_operations += 3
        
        # 1) Produce main figure for virology study.
        print(f"\n1. Generating virology research plot...")
        if plot_lovo_for_virology_research(lovo_data_for_plotting, "results/lovo_virology_analysis.png"):
            success_count += 1
            print("✓ Virology research plot generated")
        else:
            print("✗ Virology research plot failed")
        
        # 2) Produce robust LOVO box plots.
        print(f"\n2. Generating robust LOVO boxplot...")
        if plot_lovo_boxplot(lovo_data_for_plotting, "results/lovo_robust_analysis.png"):
            success_count += 1
            print("✓ Robust LOVO boxplot generated")
        else:
            print("✗ Robust LOVO boxplot failed")
        
        # 3) Produce simplified performance plots (backup).
        print(f"\n3. Generating simplified performance plot...")
        try:
            plot_simple_lovo_summary(lovo_data_for_plotting, "results/lovo_simple_summary.png")
            success_count += 1
            print("✓ Simple summary plot generated")
        except Exception as e:
            print(f"✗ Simple summary plot failed: {e}")
        
        # 4) Produce performance summary table.
        print(f"\n3. Generating performance summary...")
        if generate_lovo_summary_table(lovo_data_for_plotting, "results/lovo_performance_summary.csv"):
            success_count += 1
            print("✓ Performance summary generated")
        else:
            print("✗ Performance summary failed")
        
        # 5) Key findings.
        try:
            print_lovo_key_findings(lovo_data_for_plotting)
            print("✓ Key findings summary generated")
        except Exception as e:
            print(f"✗ Key findings failed: {e}")
        
    else:
        print(" No LOVO results available for visualization")

    # Final statistics.
    print(f"\n{'='*50}")
    print(f"PIPELINE COMPLETED: {success_count}/{total_operations} operations successful")
    if success_count == total_operations:
        print(" All operations completed successfully!")
    elif success_count > total_operations * 0.8:
        print(" Most operations completed successfully")
    else:
        print(" Some operations failed - check logs above")
    print(f"{'='*50}")

    # --- 3) UMAP projection ---
    print("\n--- UMAP Dimensionality Reduction ---")
    
    umap_configs = [
        ('3mer_', "UMAP of 3-mer Features", "results/umap_3mer_features.png"),
        ('6mer_', "UMAP of 6-mer Features", "results/umap_6mer_features.png"),
        ('RSCU_', "UMAP of RSCU Features", "results/umap_rscu_features.png"),
        (('CTDC', 'CTDT', 'CTDD', 'CTriad', 'DistancePair', 'PseAAC'), 
         "UMAP of Protein Features", "results/umap_protein_features.png"),
        (('3mer_', 'RSCU_', 'CTriad'), 
         "UMAP of Mixed Features (3-mer + RSCU + CTriad)", "results/umap_mixed_features.png")
    ]
    
    for prefix, title, output_file in umap_configs:
        total_operations += 1
        if isinstance(prefix, tuple):
            # Multiple-prefix case.
            cols = [c for c in df_for_viz.columns if any(c.startswith(p) for p in prefix)]
        else:
            # Single-prefix case.
            cols = [c for c in df_for_viz.columns if c.startswith(prefix)]
        
        if len(cols) > 0:
            if plot_umap(df_for_viz, cols, title, output_file):
                success_count += 1
        else:
            print(f" No columns found for {title}")
    
    print("✓ UMAP plots completed")
    
    # --- 4) Feature-importance analysis ---
    print("\n--- Feature Importance Analysis ---")
    total_operations += 1
    if plot_feature_importance_comparison(df_for_viz, feature_sets, "results/feature_importance"):
        success_count += 1
        print("✓ Feature importance plots completed")
    else:
        print("✗ Feature importance plots failed")
    
    # Summary.
    print(f"\n--- Visualization Summary ---")
    print(f"Total operations: {total_operations}")
    print(f"Successful operations: {success_count}")
    print(f"Success rate: {success_count/total_operations*100:.1f}%" if total_operations > 0 else "N/A")
    
    # --- 3) k-mer/RSCU clustermaps ----
    print("\n--- Clustering Heatmaps ---")
   
    all_species_order = sorted(df_for_viz['species_std'].unique()) # Ensure consistent species order.

    # 1) k-mer clustermap.
    for k in range(2, 7):  # Adjust k-mer range here (e.g., range(2,7) or range(3,5)).
        kmer_cols = [col for col in df_for_viz.columns if col.startswith(f"{k}mer_")]
        if not kmer_cols:
            print(f"[Warning] No {k}-mer features found")
            continue
        print(f"Processing {k}-mer heatmaps...")
        # Compute features for human/nonhuman separately; find shared valid features.
        human_data = df_for_viz[df_for_viz['final_label'] == 'human']
        nonhuman_data = df_for_viz[df_for_viz['final_label'] == 'nonhuman']
        
        if human_data.shape[0] == 0 or nonhuman_data.shape[0] == 0:
            print(f"[Warning] Missing data for human or nonhuman in {k}-mer")
            continue
            
        # Compute group means.
        human_means = human_data.groupby('species_std')[kmer_cols].mean()
        human_means = human_means.reindex(all_species_order)
        nonhuman_means = nonhuman_data.groupby('species_std')[kmer_cols].mean()
        nonhuman_means = nonhuman_means.reindex(all_species_order)
        # Find features valid in both groups (no NaN/Inf).
        human_clean = human_means.replace([np.inf, -np.inf], np.nan)
        nonhuman_clean = nonhuman_means.replace([np.inf, -np.inf], np.nan)
        # Find columns that are not all-NaN in either group.
        human_valid_cols = human_clean.columns[~human_clean.isnull().all()]
        nonhuman_valid_cols = nonhuman_clean.columns[~nonhuman_clean.isnull().all()]
        common_valid_cols = list(set(human_valid_cols) & set(nonhuman_valid_cols))

        if len(common_valid_cols) == 0:
            print(f"[Warning] No common valid features for {k}-mer")
            continue
    
        # Use common valid features to compute variance on combined data.
        combined_means = pd.concat([
            human_means[common_valid_cols], 
            nonhuman_means[common_valid_cols]
        ]).dropna(axis=0, how='any').dropna(axis=1, how='any')

        # Select features with highest variance.
        if combined_means.shape[1] > 0:
            final_features = top_var_features(combined_means, topn=30)
            final_cols = final_features.columns.tolist()
        
            # Cluster final features to determine order.
            if len(final_cols) > 1:
                # Standardise data.
                col_means = final_features.mean(axis=0)
                col_stds = final_features.std(axis=0).replace(0, 1)
                final_features_norm = final_features.sub(col_means, axis=1).div(col_stds, axis=1)
                # Get clustered column order.
                temp_g = sns.clustermap(final_features_norm, row_cluster=False, col_cluster=True)
                col_order = temp_g.dendrogram_col.reordered_ind
                plt.close('all')  # Ensure temporary figures are closed.
            
                # Reorder columns by clustered order.
                ordered_cols = [final_cols[i] for i in col_order]
            else:
                ordered_cols = final_cols
        else:
            print(f"[Warning] No valid features after cleaning for {k}-mer")
            continue
            
        # Plot a heatmap per label using identical columns.
        for label in ['human', 'nonhuman']:
            if label == 'human':
                data_to_plot = human_means[ordered_cols]
            else:
                data_to_plot = nonhuman_means[ordered_cols]
        
            print(f"[Info] {k}-mer {label}: Using {len(ordered_cols)} features")
            plot_clustermap(
                data_to_plot,
                f"Mean {k}-mer Features per Virus Species ({label})",
                f"mean_{k}mer_clustermap_top30_{label}.png",
                cmap='magma',
                col_cluster=False  # Disable column clustering; use the predefined order.
            )

    # 2) RSCU features — same procedure.
    rscu_cols = [col for col in df_for_viz.columns if col.startswith("RSCU_")]

    if rscu_cols:
        # Compute features for human/nonhuman separately; find shared valid features. 
        human_data = df_for_viz[df_for_viz['final_label'] == 'human']
        nonhuman_data = df_for_viz[df_for_viz['final_label'] == 'nonhuman']
        
        if human_data.shape[0] > 0 and nonhuman_data.shape[0] > 0:
            # Compute group means.
            human_means = human_data.groupby('species_std')[rscu_cols].mean()
            human_means = human_means.reindex(all_species_order)
            nonhuman_means = nonhuman_data.groupby('species_std')[rscu_cols].mean()
            nonhuman_means = nonhuman_means.reindex(all_species_order)
        
            # Find features valid in both groups (no NaN/Inf).
            human_clean = human_means.replace([np.inf, -np.inf], np.nan)
            nonhuman_clean = nonhuman_means.replace([np.inf, -np.inf], np.nan)
        
            # Find columns that are not all-NaN in either group.
            human_valid_cols = human_clean.columns[~human_clean.isnull().all()]
            nonhuman_valid_cols = nonhuman_clean.columns[~nonhuman_clean.isnull().all()]
            common_valid_cols = list(set(human_valid_cols) & set(nonhuman_valid_cols))
        
            if len(common_valid_cols) > 0:
                # Use common valid features to compute variance on combined data.
                combined_means = pd.concat([
                    human_means[common_valid_cols], 
                    nonhuman_means[common_valid_cols]
                ]).dropna(axis=0, how='any').dropna(axis=1, how='any')
                
                # Select features with highest variance.
                if combined_means.shape[1] > 0:
                    final_features = top_var_features(combined_means, topn=30)
                    final_cols = final_features.columns.tolist()
                    
                    # Cluster final features to determine order.
                    if len(final_cols) > 1:
                        # Standardise data.
                        col_means = final_features.mean(axis=0)
                        col_stds = final_features.std(axis=0).replace(0, 1)
                        final_features_norm = final_features.sub(col_means, axis=1).div(col_stds, axis=1)
                        
                        # Get clustered column order.
                        temp_g = sns.clustermap(final_features_norm, row_cluster=False, col_cluster=True)
                        col_order = temp_g.dendrogram_col.reordered_ind
                        plt.close('all')  # Ensure temporary figures are closed.
                    
                        # Reorder columns by clustered order.
                        ordered_cols = [final_cols[i] for i in col_order]
                    else:
                        ordered_cols = final_cols
                
                    for label in ['human', 'nonhuman']:
                        if label == 'human':
                            data_to_plot = human_means[ordered_cols]
                        else:
                            data_to_plot = nonhuman_means[ordered_cols]    
                            
                        print(f"[Info] RSCU {label}: Using {len(ordered_cols)} features")
                        plot_clustermap(
                            data_to_plot,
                            f"Mean RSCU Features per Virus Species ({label})",
                            f"mean_rscu_clustermap_top30_{label}.png",
                            cmap='viridis',
                            col_cluster=False  
                        )
                else:
                    print("[Warning] No valid features after cleaning for RSCU")
            else:
                print("[Warning] No common valid features for RSCU")
        else:
            print("[Warning] Missing data for human or nonhuman in RSCU")
    else:
        print("[Warning] No RSCU features found")

    print("All heatmaps finished.")

def _pick_first_existing(paths):
    """Given a list of candidate paths, return the first existing absolute path;
    if none exists, return the absolute path of the first candidate (used in error messages)."""
    for p in paths:
        if p and Path(p).exists():
            return str(Path(p).resolve())
    return str(Path(paths[0]).resolve())

# Main entry.
if __name__ == '__main__':
    # Determine repository root based on script location (one level above src)
    SRC_DIR  = Path(__file__).resolve().parent
    REPO_DIR = SRC_DIR.parent
    DATA_DIR = REPO_DIR / "data"
    RESULTS_DIR = REPO_DIR / "results"

    parser = argparse.ArgumentParser()
    # Note: --nucleotide / --nonhuman correspond to "human / nonhuman nucleotide FASTA"
    parser.add_argument("--nucleotide", help="Human nucleotide FASTA (legacy: sequences_human.fasta; alt: data/nucleotide.fasta)")
    parser.add_argument("--nonhuman",  help="Nonhuman nucleotide FASTA (legacy: sequences_nothuman.fasta; alt: data/nonhuman.fasta)")
    parser.add_argument("--cds",       help="CDS FASTA (legacy: sequences_CDS.fasta; alt: data/cds.fasta)")
    parser.add_argument("--protein",   help="Protein FASTA (legacy: sequences_protein.fasta; alt: data/protein.fasta)")
    parser.add_argument("--reference", help="CAI reference FASTA (legacy: human_HK_CDS.cleaned.fasta; alt: data/human_HK_CDS.cleaned.fasta)")
    parser.add_argument("--outdir",    default=str(RESULTS_DIR), help="Output directory (default: <repo>/results)")
    args = parser.parse_args()

    # Auto-detection: prioritize legacy naming (repo root), fallback to README naming (data/)
    human_fasta = _pick_first_existing([
        args.nucleotide or (REPO_DIR / "sequences_human.fasta"),
        DATA_DIR / "nucleotide.fasta"
    ])
    nonhuman_fasta = _pick_first_existing([
        args.nonhuman or (REPO_DIR / "sequences_nothuman.fasta"),
        DATA_DIR / "nonhuman.fasta"
    ])
    cds_fasta = _pick_first_existing([
        args.cds or (REPO_DIR / "sequences_CDS.fasta"),
        DATA_DIR / "cds.fasta"
    ])
    protein_fasta = _pick_first_existing([
        args.protein or (REPO_DIR / "sequences_protein.fasta"),
        DATA_DIR / "protein.fasta"
    ])
    reference_fasta = _pick_first_existing([
        args.reference or (REPO_DIR / "human_HK_CDS.cleaned.fasta"),
        DATA_DIR / "human_HK_CDS.cleaned.fasta"
    ])

    # ---- Record a snapshot of root directory files *before* execution, for collecting scattered outputs later ----
    before_files = {p.resolve() for p in REPO_DIR.glob("*") if p.is_file()}

    # Build absolute-path config
    OUTDIR = Path(args.outdir).resolve()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    config = {
        'human_fasta': human_fasta,
        'nonhuman_fasta': nonhuman_fasta,
        'cds_fasta': cds_fasta,
        'protein_fasta': protein_fasta,
        'reference_fasta': reference_fasta,
        'outdir': str(OUTDIR)
    }

    log_path = OUTDIR / 'pipeline_log.txt'
    log_file = open(log_path, 'w', encoding='utf-8')
    sys.stdout = Tee(sys.__stdout__, log_file)
    sys.stderr = Tee(sys.__stderr__, log_file)

    print("=== Script Start ===")
    print("[INFO] Repository root:", REPO_DIR)
    print("[INFO] Using config:")
    for k, v in config.items():
        print(f"  - {k}: {v}")

    # Run main pipeline (keep CWD at repo root to ensure reference files with relative paths can be read)
    run_pipeline(config)

    # ---- After execution: collect "newly generated files scattered in repo root" into results ----
    after_files = {p.resolve() for p in REPO_DIR.glob("*") if p.is_file()}
    new_files = after_files - before_files

    # Whitelist of input filenames (prevent them from being moved)
    input_basenames = {
        'sequences_human.fasta',
        'sequences_nothuman.fasta',
        'sequences_CDS.fasta',
        'sequences_protein.fasta',
        'human_HK_CDS.cleaned.fasta',
        'Orthobunyavirus_taxonomy.xlsx',
        'ICTV_2024_Taxa Renamed or Abolished.csv',
        'ICTV_2024_MSL.csv',
        'ICTV Orthobunyavirus.csv',
        'human_virus_DB_species_clean.txt',
    }

    # Only relocate "newly added top-level files", restricted to common result suffixes, excluding input lists and log files just written
    allowed_suffixes = {'.csv', '.tsv', '.txt', '.png', '.pdf', '.svg', '.json', '.log'}
    moved = []
    for f in new_files:
        if (f.parent == REPO_DIR
            and f.name not in input_basenames
            and f.suffix.lower() in allowed_suffixes
            and f != log_path.resolve()):
            dest = OUTDIR / f.name
            try:
                shutil.move(str(f), str(dest))
                moved.append(dest.name)
            except Exception as e:
                print(f"[WARN] Failed to move {f.name} to results: {e}")

    if moved:
        print(f"[INFO] Top-level outputs moved to {OUTDIR.name}/: {', '.join(moved)}")
    else:
        print("[INFO] No new top-level outputs to move.")

    print("=== Script End ===")
    log_file.close()
