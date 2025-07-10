#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import shlex
import tempfile
import shutil
import json
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
import Bio
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


############################################################
#                   RIPPredictor Class
############################################################

class RIPPredictor:
    """
    使用Bakta对输入fasta进行注释，然后使用MMseqs2和HMMER搜索标记RIP蛋白。
    """
    def __init__(self, fasta_file, db_file, hmm_file, output_dir, e_value, min_score, min_alnlen):
        self.fasta_file = fasta_file
        self.db_file = db_file
        self.hmm_file = hmm_file
        self.e_value = e_value
        self.min_score = min_score
        self.min_alnlen = min_alnlen
        self.output_dir = output_dir
        self.gbff_file = None

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)
        self.output_folder = os.path.join(self.output_dir, os.path.splitext(os.path.basename(self.fasta_file))[0])
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)

    def predict(self):
        """
        运行RIP预测流程
        Returns:
            gbff_file (str): 生成的gbff文件路径
            results_df (pd.DataFrame): 完整注释结果DataFrame
        """
        ### START: lines that I edited ###
        # self._run_bakta()
        self.gbff_file = f"../bakta-annotations/results/{os.path.splitext(os.path.basename(self.fasta_file))[0]}/{os.path.splitext(os.path.basename(self.fasta_file))[0]}.gbff"
        ### END ###
        genes_df = self._get_gene_info()
        blast_results = self._run_mmseqs_search(genes_df)
        hmm_results = self._run_hmm_search(genes_df)
        results_df = self._process_results(genes_df, blast_results, hmm_results)
        self._update_gbff(results_df)
        return self.gbff_file, results_df

    def _run_bakta(self):
        """
        使用bakta对输入fasta进行基因组注释
        """
        output_folder = os.path.join(self.output_dir, os.path.splitext(os.path.basename(self.fasta_file))[0])
        command = [
            "bakta",
            "--complete",
            "--db",
            "Bakta/db-light",
            "--skip-plot",
            "--keep-contig-headers",
            "--output",
            output_folder,
            self.fasta_file
        ]
        subprocess.run(command, check=True)
        self.gbff_file = os.path.join(
            output_folder,
            os.path.splitext(os.path.basename(self.fasta_file))[0] + ".gbff"
        )
        if not os.path.exists(self.gbff_file):
            raise FileNotFoundError("gbff file not found after running bakta.")

    def _get_gene_info(self):
        """
        从gbff文件中提取基因信息
        Returns:
            genes_df (pd.DataFrame): 基因信息的DataFrame
        """
        genes_df = pd.DataFrame(columns=[
            'Accession_Number', 'Sequence_length', 'gene_start', 'gene_end',
            'gene_strand', 'translation', 'gene_order', 'gene_id', 'gene', 'product'
        ])

        line = 0
        for record in SeqIO.parse(self.gbff_file, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    if isinstance(feature.location, CompoundLocation):
                        if feature.strand == -1:
                            gene_start = int(feature.location.parts[-1].start)
                            gene_end = int(feature.location.parts[0].end)
                        else:
                            gene_start = int(feature.location.parts[0].start)
                            gene_end = int(feature.location.parts[-1].end)
                    else:
                        gene_start = int(feature.location.start)
                        gene_end = int(feature.location.end)

                    genes_df.at[line, 'Accession_Number'] = record.id
                    genes_df.at[line, 'Sequence_length'] = len(record.seq)
                    genes_df.at[line, 'gene_strand'] = str(feature.strand)
                    genes_df.at[line, 'translation'] = str(feature.qualifiers.get('translation', [''])[0])
                    genes_df.at[line, 'gene_order'] = line + 1
                    genes_df.at[line, 'gene_id'] = feature.qualifiers.get('locus_tag', [''])[0]
                    genes_df.at[line, 'gene'] = feature.qualifiers.get('gene', ['-'])[0].lower()
                    genes_df.at[line, 'product'] = str(feature.qualifiers.get('product', [''])[0])
                    genes_df.at[line, 'gene_start'] = gene_start
                    genes_df.at[line, 'gene_end'] = gene_end
                    line += 1
        return genes_df

    def _run_mmseqs_search(self, genes_df):
        """
        使用MMseqs2对提取的基因序列进行搜索
        """
        fasta_folder = os.path.join(self.output_dir, "fasta_files")
        os.makedirs(fasta_folder, exist_ok=True)
        all_cds_fasta = os.path.join(fasta_folder, "all_cds.fasta")

        with open(all_cds_fasta, 'w') as fasta_file:
            for _, row in genes_df.iterrows():
                if row['translation']:
                    fasta_file.write(f">{row['gene_id']}_{row['gene_order']}\n{row['translation']}\n")

        blast_results = []
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "aln.tsv")

            cmd = f"mmseqs easy-search {shlex.quote(all_cds_fasta)} {shlex.quote(self.db_file)} " \
                  f"{shlex.quote(output_file)} {shlex.quote(os.path.join(temp_dir, 'tmp'))} " \
                  f"--format-output 'query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits'"
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

            if os.path.exists(output_file):
                results_dict = {}
                with open(output_file, "r") as results:
                    for line in results:
                        fields = line.strip().split("\t")
                        if len(fields) == 12:
                            query, target, pident, alnlen, qstart, qend, qlen, tstart, tend, tlen, evalue, bits = fields
                            gene_id = query.rsplit('_', 1)[0]
                            score = float(bits)
                            current_evalue = float(evalue)
                            qlen = int(qlen)
                            tlen = int(tlen)
                            alnlen = int(alnlen)
                            if (current_evalue <= self.e_value and
                                score >= self.min_score and
                                qlen >= 140 and
                                0.5 * tlen <= qlen <= 1.5 * tlen and
                                alnlen >= self.min_alnlen):
                                if gene_id not in results_dict or score > results_dict[gene_id][6]:
                                    results_dict[gene_id] = [
                                        gene_id, target, float(pident), alnlen, qlen,
                                        tlen, current_evalue, score
                                    ]
                blast_results = list(results_dict.values())

        return pd.DataFrame(blast_results, columns=[
            "Gene_ID", "mmseqs_hit", "mmseqs_Identity",
            "mmseqs_Alignmentlength", "mmseqs_Querylength",
            "mmseqs_Subjectlength", "mmseqs_Evalue", "mmseqs_Bitscore"
        ])

    def _parse_hmm_file(self):
        hmm_info = {}
        with open(self.hmm_file, "r") as f:
            current_hmm = None
            for line in f:
                if line.startswith("NAME"):
                    current_hmm = line.split()[1]
                    hmm_info[current_hmm] = {"length": 0, "tc": None, "ga": None, "nc": None}
                elif line.startswith("LENG") and current_hmm:
                    hmm_info[current_hmm]["length"] = int(line.split()[1])
                elif line.startswith("TC") and current_hmm:
                    hmm_info[current_hmm]["tc"] = float(line.split()[1])
                elif line.startswith("GA") and current_hmm:
                    hmm_info[current_hmm]["ga"] = float(line.split()[1])
                elif line.startswith("NC") and current_hmm:
                    hmm_info[current_hmm]["nc"] = float(line.split()[1])
        return hmm_info

    def _run_hmm_search(self, genes_df):
        """
        使用HMMER对基因序列进行搜索
        """
        hmm_results = []
        hmm_info = self._parse_hmm_file()

        with tempfile.TemporaryDirectory() as temp_dir:
            for _, row in genes_df.iterrows():
                if row['translation']:
                    sequence_id = f"{row['gene_id']}_{row['gene_order']}"
                    result = self._run_hmmscan(row['translation'], sequence_id, row['gene_id'], temp_dir, hmm_info)
                    if result:
                        hmm_results.append(result)

        return pd.DataFrame(hmm_results, columns=[
            "gene_id", "hmm_hit", "hmm_evalue", "hmm_alignment", "hmm_score", "hmm_accuracy"
        ])

    def _run_hmmscan(self, sequence, sequence_id, gb_id, temp_dir, hmm_info):
        safe_id = "".join(c for c in sequence_id if c.isalnum() or c in "._- ")
        temp_fasta = os.path.join(temp_dir, f"{safe_id}.fasta")
        temp_domtblout = os.path.join(temp_dir, f"{safe_id}.domtblout")

        with open(temp_fasta, 'w') as f:
            f.write(f">{sequence_id}\n{sequence}\n")

        cmd = f"hmmscan --noali --notextw --domtblout {shlex.quote(temp_domtblout)} {shlex.quote(self.hmm_file)} {shlex.quote(temp_fasta)}"
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        best_result = None
        best_score = float('-inf')
        with open(temp_domtblout, "r") as domtblout:
            for line in domtblout:
                if not line.startswith("#"):
                    fields = line.strip().split()
                    if len(fields) >= 23:
                        hmm_name = fields[0]
                        e_value = float(fields[6])
                        score = float(fields[13])
                        ali_from = int(fields[17])
                        ali_to = int(fields[18])
                        acc = float(fields[21])
                        score_threshold = hmm_info.get(hmm_name, {}).get("nc", 0)
                        if score >= score_threshold and score > best_score:
                            best_score = score
                            alignment = f"{ali_from}-{ali_to}"
                            best_result = (gb_id, hmm_name, e_value, alignment, score, acc)
        return best_result

    def _process_results(self, genes_df, blast_results, hmm_results):
        """
        合并MMseqs2和HMMER结果，并标记RIP蛋白
        """
        combined_df = pd.merge(
            genes_df,
            blast_results,
            how='left',
            left_on='gene_id',
            right_on='Gene_ID'
        )
        combined_df = pd.merge(
            combined_df,
            hmm_results,
            how='left',
            on='gene_id'
        )
        combined_df['RIP'] = combined_df.apply(
            lambda row: 1 if pd.notna(row['mmseqs_hit']) or pd.notna(row['hmm_hit']) else 0, axis=1
        )

        # bakta_csv_path = os.path.join(self.output_dir, "bakta.csv")
        # combined_df.to_csv(bakta_csv_path, index=False)

        rip_df = combined_df[combined_df['RIP'] == 1].copy()
        hmm_cols = ["hmm_hit", "hmm_evalue", "hmm_alignment", "hmm_score", "hmm_accuracy"]
        rip_df.drop(columns=hmm_cols, inplace=True, errors='ignore')
        rip_csv_path = os.path.join(self.output_dir, "RIP.csv")
        rip_df.to_csv(rip_csv_path, index=False)

        return combined_df

    def _update_gbff(self, results_df):
        """
        在gbff文件中更新RIP蛋白标记
        """
        records = []
        for record in SeqIO.parse(self.gbff_file, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS' and 'locus_tag' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    gene_row = results_df[results_df['gene_id'] == locus_tag]
                    if not gene_row.empty and gene_row['RIP'].iloc[0] == 1:
                        feature.qualifiers['product'] = ['RIP']
            records.append(record)

        with open(self.gbff_file, 'w') as output_handle:
            SeqIO.write(records, output_handle, 'genbank')


############################################################
#                 Iteron Class
############################################################

class Iteron:
    """
    对intergenic序列进行iteron分析。
    """
    def __init__(self, intergenic_sequences, min_iteron, max_iteron):
        self.intergenic_sequences = intergenic_sequences
        self.min_iteron = min_iteron
        self.max_iteron = max_iteron

    def dna_sequence_score(self, seq1, seq2):
        match_score = 2
        partial_match_score = 0
        no_match_score = 0
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be equal length.")
        score = 0
        for base1, base2 in zip(seq1, seq2):
            if base1 == base2:
                score += match_score
            elif (base1 == 'A' and base2 == 'T') or (base1 == 'T' and base2 == 'A') or \
                 (base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G'):
                score += partial_match_score
            else:
                score += no_match_score
        return score

    def calculate_entropy(self, sequences):
        aligned_positions = list(zip(*sequences))
        entropy_scores = []
        sequence_length = len(aligned_positions)
        for position in aligned_positions:
            counts = np.unique(position, return_counts=True)
            total = sum(counts[1])
            probabilities = [count / total for count in counts[1] if '-' not in counts[0]]
            entropy = -sum(p * np.log2(p) for p in probabilities)
            entropy_scores.append(entropy)
        ebcons = max(0, (sum(1 - e for e in entropy_scores) / sequence_length) * 100)
        return entropy_scores, ebcons

    def clustering_index_std(self, distances):
        if len(distances) <= 1:
            return 1
        total_span = sum(distances)
        clustering_index = 1 / (1 + np.log(1 + total_span/len(distances)))
        return clustering_index

    def process_sequence(self, sequence):
        downstream_range = 100
        all_kmer_results = []
        processed_regions = set()

        for kmer in range(self.max_iteron, self.min_iteron, -1):
            min_score = int(kmer * 0.75 * 2)
            sequence_length = len(sequence)
            results = []

            for start in range(sequence_length - kmer + 1):
                if any(start <= r[0] <= start + kmer or start <= r[1] <= start + kmer for r in processed_regions):
                    continue
                kmer_seq = sequence[start:start + kmer]
                matching_sequences_positions = []
                current_start = start + kmer

                while current_start < sequence_length:
                    downstream_seq = sequence[current_start:current_start + downstream_range]
                    best_match = None
                    best_score = 0
                    best_distance = float('inf')

                    for downstream_start in range(len(downstream_seq) - kmer + 1):
                        candidate_seq = downstream_seq[downstream_start:downstream_start + kmer]
                        score = self.dna_sequence_score(kmer_seq, candidate_seq)
                        if score >= min_score:
                            if score > best_score or (score == best_score and downstream_start < best_distance):
                                best_match = (candidate_seq, current_start + downstream_start, score, downstream_start)
                                best_score = score
                                best_distance = downstream_start

                    if best_match is not None:
                        matching_sequences_positions.append(best_match)
                        current_start = best_match[1] + kmer
                    else:
                        break

                if len(matching_sequences_positions) >= 2:
                    total_score = sum(m[2] for m in matching_sequences_positions)
                    normalized_score = total_score / (len(matching_sequences_positions) * kmer * 2)
                    distances = [matching_sequences_positions[0][1] - (start + kmer)]
                    distances.extend([
                        matching_sequences_positions[i+1][1] - (matching_sequences_positions[i][1] + kmer)
                        for i in range(len(matching_sequences_positions) - 1)
                    ])
                    clustering_index = self.clustering_index_std(distances)
                    aligned_sequences = [kmer_seq] + [m[0] for m in matching_sequences_positions]
                    entropy_scores, ebcons = self.calculate_entropy(aligned_sequences)

                    result_dict = {
                        "kmer_length": kmer,
                        "repeat_count": len(matching_sequences_positions) + 1,
                        "total_position": f"{start}-{matching_sequences_positions[-1][1] + kmer}",
                        "match_score": normalized_score,
                        "score": (len(matching_sequences_positions) + 1) * normalized_score,
                        "kmer": kmer_seq,
                        "clustering_index": clustering_index,
                        "entropy_score": ebcons,
                        "positions": [
                            {"match_sequence": kmer_seq, "position_range": f"{start}-{start + kmer}"}
                        ] + [
                            {"match_sequence": mm[0], "position_range": f"{mm[1]}-{mm[1] + kmer}"}
                            for mm in matching_sequences_positions
                        ]
                    }
                    results.append(result_dict)
                    processed_regions.add((start, matching_sequences_positions[-1][1] + kmer))
            all_kmer_results.extend(results)
        return all_kmer_results

    def filter_kmer_results(self, all_kmer_results):
        final_filtered = []
        for result in all_kmer_results:
            final_score = result['kmer_length'] * result['clustering_index'] * \
                          result['repeat_count'] * (result['entropy_score'] / 100)
            result['final_score'] = final_score
            final_filtered.append(result)
        final_filtered.sort(key=lambda x: x['final_score'], reverse=True)
        return final_filtered[0] if final_filtered else None

    def analyze(self, sequence):
        processed_results = self.process_sequence(sequence)
        filtered_result = self.filter_kmer_results(processed_results)

        if filtered_result:
            iteron_score = filtered_result.get('final_score', 0)
            iteron_entropy_score = filtered_result.get('entropy_score', 0)
            positions = filtered_result.get('positions', [])
            position_data = []
            for p in positions:
                s, e = p["position_range"].split("-")
                position_data.append({"start": int(s), "end": int(e)})
            iteron_df = pd.DataFrame(position_data)
            return iteron_score, iteron_df, iteron_entropy_score
        else:
            return 0, pd.DataFrame(columns=["start", "end"]), 0


############################################################
#                ZCurveAnalyzer Class
############################################################

class ZCurveAnalyzer:
    """
    对intergenic序列进行ZCurve分析，从中计算AT score。
    """
    def __init__(self, intergenic_sequences, sigma):
        self.intergenic_sequences = intergenic_sequences
        self.sigma = sigma

    @staticmethod
    def regression_line(x_array, y_array):
        x_array = np.array(x_array)
        y_array = np.array(y_array)
        x_array_mean = np.mean(x_array)
        y_array_mean = np.mean(y_array)
        Sxx = sum((xi - x_array_mean)**2 for xi in x_array)
        Sxy = sum((x_array[i] - x_array_mean)*(y_array[i] - y_array_mean) for i in range(len(x_array)))
        k = Sxy / Sxx
        return k

    def Z_prime_curves(self, seq):
        length = len(seq)
        para_array = np.zeros((3, length), dtype='float64')
        sum_a = sum_c = sum_g = sum_t = 0

        for i, base in enumerate(seq.upper()):
            para_array[0, i] = i + 1
            if base == 'A':
                sum_a += 1
            elif base == 'C':
                sum_c += 1
            elif base == 'G':
                sum_g += 1
            elif base == 'T':
                sum_t += 1
            para_array[1, i] = (sum_a + sum_t) - (sum_g + sum_c)

        X = para_array[0]
        Y = para_array[1]
        k_para = self.regression_line(X, Y)

        for n_para in para_array[0]:
            index = int(n_para - 1)
            para_array[2, index] = para_array[1, index] - float(k_para) * n_para

        df_para = pd.DataFrame(np.transpose(para_array),
                               columns=['n', 'Zn_curve', 'Zn_prime_curve'])
        df_para = df_para.astype({"n": int})
        df_para['Zn_prime_curve_smooth'] = gaussian_filter1d(df_para['Zn_prime_curve'], self.sigma)
        return df_para

    def analyze(self, sequence):
        df = self.Z_prime_curves(sequence)
        peak_indices = argrelextrema(df['Zn_prime_curve_smooth'].values, np.greater)[0]
        valley_indices = argrelextrema(-df['Zn_prime_curve_smooth'].values, np.greater)[0]
        split_indices = sorted(np.concatenate((peak_indices, valley_indices)))
        split_indices = np.concatenate(([0], split_indices, [len(sequence) - 1]))

        max_at_score = -float('inf')
        max_segment = None
        for i in range(len(split_indices) - 1):
            start_n = df.loc[split_indices[i], 'n']
            end_n = df.loc[split_indices[i + 1], 'n']
            subseq_length = end_n - start_n - 1
            if subseq_length <= 0:
                continue
            subseq = sequence[start_n + 1:end_n - 1]
            at_content = (subseq.count('A') + subseq.count('T')) / subseq_length
            x = df[(df['n'] >= start_n) & (df['n'] <= end_n)]['n']
            y = df[(df['n'] >= start_n) & (df['n'] <= end_n)]['Zn_prime_curve_smooth']
            slope = self.regression_line(x, y)
            at_score = subseq_length * slope * at_content
            if at_score > max_at_score:
                max_at_score = at_score
                max_segment = (start_n, end_n, subseq_length, slope, at_content, at_score)

        if max_segment:
            at_df = pd.DataFrame([{
                'start_n': max_segment[0],
                'end_n': max_segment[1],
                'slope': max_segment[3],
                'at_content': max_segment[4],
                'at_score': max_segment[5]
            }])
            return max_at_score, at_df
        else:
            at_df = pd.DataFrame(columns=['start_n', 'end_n', 'slope', 'at_content', 'at_score'])
            return 0, at_df


class RNAMotifAnalyzer:
    """
    使用cmscan对RNA Motif进行分析
    """
    def __init__(self, cm_folder="../orivfinder/app/RNA_db", e_value_threshold=1e-5):
        self.cm_folder = cm_folder
        self.e_value_threshold = e_value_threshold
        self.temp_dir = tempfile.mkdtemp()
        self.cm_files = [f for f in os.listdir(self.cm_folder) if f.endswith('.cm')]
        # 明确指定哪些文件使用--cut_tc，哪些使用E-value
        self.cut_tc_files = ["PAGEV.cm", "RNAI.cm"]
        self.evalue_files = ["RNAII.cm"]
        

    def _create_temp_fasta(self, sequence):
        temp_fasta = os.path.join(self.temp_dir, "temp.fasta")
        with open(temp_fasta, 'w') as f:
            f.write(">temp_seq\n")
            f.write(sequence + "\n")
        return temp_fasta

    def _run_cmscan(self, fasta_file, cm_file):
        output_file = os.path.join(self.temp_dir, f"cmscan_output_{os.path.splitext(os.path.basename(cm_file))[0]}.txt")
        
        # 根据文件名决定使用哪种参数
        if cm_file in self.cut_tc_files:

            cmd = [
                'cmscan',
                '--noali',
                '-g',
                '--cut_tc',
                '--nohmmonly',
                '--cpu', '1',
                '--tblout', output_file,
                os.path.join(self.cm_folder, cm_file),
                fasta_file
            ]
        else:  # 默认使用E-value或明确指定的E-value文件

            cmd = [
                'cmscan',
                '--noali',
                '-g',
                '-E', str(self.e_value_threshold),
                '--nohmmonly',
                '--cpu', '1',
                '--tblout', output_file,
                os.path.join(self.cm_folder, cm_file),
                fasta_file
            ]
        
        try:

            subprocess.run(cmd, check=True, capture_output=True, text=True)

        except subprocess.CalledProcessError as e:

            # 如果出错，创建一个空的输出文件以便后续处理
            with open(output_file, 'w') as f:
                f.write("# No results due to error\n")
        
        return output_file

    def _parse_cmscan_output(self, output_file):
        results = []
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 16:  # 确保有足够的列
                            start = int(parts[7])
                            end = int(parts[8])
                            if start > end:
                                start, end = end, start
                            results.append({
                                'name': parts[0],
                                'accession': parts[1],
                                'start': start,
                                'end': end,
                                'strand': parts[9],
                                'evalue': float(parts[15])
                            })
        except Exception as e:
            print(f" {output_file} erro: {str(e)}")
        
        return results

    def _filter_overlapping_hits(self, hits):
        sorted_hits = sorted(hits, key=lambda x: x['evalue'])
        filtered_hits = []
        covered_regions = []
        for hit in sorted_hits:
            start = min(hit['start'], hit['end'])
            end = max(hit['start'], hit['end'])
            current_range = (start, end)
            is_overlapping = False
            for region in covered_regions:
                if not (current_range[1] < region[0] or current_range[0] > region[1]):
                    is_overlapping = True
                    break
            if not is_overlapping:
                filtered_hits.append(hit)
                covered_regions.append(current_range)
        return filtered_hits

    def analyze(self, sequence):
        all_results = []
        temp_fasta = self._create_temp_fasta(sequence)

        try:

            for i, cm_file in enumerate(self.cm_files):

                output_file = self._run_cmscan(temp_fasta, cm_file)
                results = self._parse_cmscan_output(output_file)

                all_results.extend(results)
        except Exception as e:
            print(f"{str(e)}")
        finally:
            if os.path.exists(temp_fasta):
                os.remove(temp_fasta)


        print(f"{len(all_results)} match")
        filtered_results = self._filter_overlapping_hits(all_results)
        print(f" {len(filtered_results)} ")

        if filtered_results:
            RNA_df = pd.DataFrame(filtered_results)[['name', 'start', 'end', 'strand', 'evalue']]
            RNA_df['evalue'] = RNA_df['evalue'].apply(lambda x: f"{x:.2e}")
            RNA_score = 1000
        else:
            RNA_df = pd.DataFrame(columns=['name', 'start', 'end', 'strand', 'evalue'])
            RNA_score = 0

        return RNA_df, RNA_score
        
    def __del__(self):
        if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
############################################################
#               PatternFinder Class
############################################################

class PatternFinder:
    """
    在intergenic序列中寻找特定DNA序列模式和RNA结构motif。
    """
    iupac_dict = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    iupac_complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N'
    }

    patterns = {
        "DnaA_box": {"sequence": "TTATCCACA", "max_mismatch": 1},
        "ctra": {"sequence": "TTAA.{7}TTAA", "regex": True},
        "Fis_binding_motif": {"sequence": "GNNNAWWWWWTNNNC", "iupac": True},
        "IHF_binding_motif": {"sequence": "WATCAANNNNTTR", "iupac": True},
        "cole2_3": {"sequence": "TATCAGATAACAGC", "max_mismatch": 1},
        "incQ_conserved": {"sequence": "ATCCACTAGGCGCGG", "max_mismatch": 2},
        "Nick site": {"sequence": "AAAACCGGCTACTCTAATAGCCGGTT", "max_mismatch": 2},    #nic 1
        "Nick site": {"sequences": [                                                   #nic 2
            {"sequence": "GTTTTTTCTTATCTTGATAC", "max_mismatch": 2},
            {"sequence": "TTTATACTTAAGGGATA", "max_mismatch": 2},
            {"sequence": "TTCTTATCTTGATA", "max_mismatch": 2},
        ]},
        "Nick site": {"sequences": [                                                   #nic3
            {"sequence": "TACTACGACCCCCCC", "max_mismatch": 2},
            {"sequence": "TACTACAACACCCCCCTATGT", "max_mismatch": 2},
            {"sequence": "ATGTAATGGATAAGAAAATTACTTTTAAAAT", "max_mismatch": 2},
            {"sequence": "GAGGGAATGAAATTCCCTCAT", "max_mismatch": 2},
            {"sequence": "TAAAAAAATACACGATTTTTTA", "max_mismatch": 2}
        ]},
    }

    def __init__(self, intergenic_df):
        if not isinstance(intergenic_df, pd.DataFrame):
            raise TypeError("intergenic_df must be a pandas DataFrame")
        self.intergenic_df = intergenic_df

    def count_mismatches(self, seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        return sum(1 for a, b in zip(seq1, seq2) if a != b)

    def find_pattern_with_mismatch(self, seq, pattern, max_mismatch, strand):
        positions = []
        pattern_len = len(pattern)
        for i in range(len(seq) - pattern_len + 1):
            sub_seq = seq[i:i + pattern_len]
            mismatches = self.count_mismatches(sub_seq, pattern)
            if mismatches <= max_mismatch:
                positions.append((i + 1, i + pattern_len, sub_seq, mismatches, strand))
        return positions

    def find_pattern_with_regex(self, seq, pattern, strand):
        positions = []
        for match in re.finditer(pattern, seq):
            start, end = match.start(), match.end()
            sub_seq = seq[start:end]
            positions.append((start + 1, end, sub_seq, 0, strand))
        return positions

    def iupac_to_regex(self, pattern):
        return ''.join(self.iupac_dict[base] for base in pattern)

    def iupac_reverse_complement(self, pattern):
        return ''.join(self.iupac_complement[base] for base in reversed(pattern))

    def find_patterns_in_sequence(self, sequence):
        from_pattern = RNAMotifAnalyzer()
        matches = []
        seen_positions = set()
        pattern_score = 0

        RNA_df, RNA_score = from_pattern.analyze(sequence)
        pattern_score += RNA_score

        scoring_rules = {
            'DnaA_box': {0: 20, 1: 10},
            'ctra': 5,
            'Fis_binding_motif': 5,
            'IHF_binding_motif': 5,
            'cole2_3': {0: -50, 1: -40, 2: -30, 3: -20},
            'incQ_conserved': {0: -50, 1: -40, 2: -30, 3: -20},
            'Nick site': {0: -50, 1: -40, 2: -30, 3: -20},
            # 'inc2': {0: -50, 1: -40, 2: -30, 3: -20},
            # 'inc3': {0: -50, 1: -40, 2: -30, 3: -20},
        }

        for pattern_name, details in self.patterns.items():
            if 'sequences' in details:
                # 多序列pattern处理
                for seq_details in details['sequences']:
                    pattern_seq = seq_details['sequence']
                    max_mismatch = seq_details['max_mismatch']
                    # Forward
                    for match in self.find_pattern_with_mismatch(sequence, pattern_seq, max_mismatch, '+'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name][match[3]]

                    # Reverse
                    reverse_pattern = str(Seq(pattern_seq).reverse_complement())
                    for match in self.find_pattern_with_mismatch(sequence, reverse_pattern, max_mismatch, '-'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name][match[3]]

            else:
                # 单序列pattern
                pattern_seq = details['sequence']
                if 'regex' in details and details['regex']:
                    # Regex匹配
                    # Forward
                    for match in self.find_pattern_with_regex(sequence, pattern_seq, '+'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name]

                    # Reverse
                    reverse_pattern = str(Seq(pattern_seq).reverse_complement())
                    for match in self.find_pattern_with_regex(sequence, reverse_pattern, '-'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name]

                elif 'iupac' in details and details['iupac']:
                    regex_pattern = self.iupac_to_regex(pattern_seq)
                    # Forward
                    for match in self.find_pattern_with_regex(sequence, regex_pattern, '+'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name]

                    # Reverse
                    reverse_pattern = self.iupac_reverse_complement(pattern_seq)
                    reverse_regex_pattern = self.iupac_to_regex(reverse_pattern)
                    for match in self.find_pattern_with_regex(sequence, reverse_regex_pattern, '-'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            pattern_score += scoring_rules[pattern_name]

                else:
                    # 普通序列匹配（带错配）
                    max_mismatch = details['max_mismatch']
                    # Forward
                    for match in self.find_pattern_with_mismatch(sequence, pattern_seq, max_mismatch, '+'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            if pattern_name in scoring_rules:
                                pattern_score += scoring_rules[pattern_name][match[3]]

                    # Reverse
                    reverse_pattern = str(Seq(pattern_seq).reverse_complement())
                    for match in self.find_pattern_with_mismatch(sequence, reverse_pattern, max_mismatch, '-'):
                        pos_key = (match[0], match[1])
                        if pos_key not in seen_positions:
                            matches.append((pattern_name, *match))
                            seen_positions.add(pos_key)
                            if pattern_name in scoring_rules:
                                pattern_score += scoring_rules[pattern_name][match[3]]

        pattern_df = pd.DataFrame(matches, columns=['name', 'start', 'end', 'matched_sequence', 'mismatch', 'strand'])
        pattern_df = pattern_df[['name', 'start', 'end', 'strand', 'mismatch']]

        return pattern_score, RNA_df, pattern_df


############################################################
#           AnnotationAnalyzer Class
############################################################



class AnnotationAnalyzer:
    """
    对gbff文件进行注释分析，包括基因提取、假设蛋白过滤、intergenic序列分析、合并特征等。
    """
    def __init__(self, gb_file, min_iteron, max_iteron, sigma):
        self.gb_file = gb_file
        self.gb_records = list(SeqIO.parse(self.gb_file, "genbank"))
        self.min_iteron = min_iteron
        self.max_iteron = max_iteron
        self.sigma = sigma  # Added for ZCurveAnalyzer

    def gbk2df(self):
        genes_df = pd.DataFrame(columns=[
            'Accession_Number', 'Sequence_length', 'gene_start', 'gene_end',
            'gene_strand', 'translation', 'gene_order', 'gene_id', 'gene', 'product'
        ], dtype=str)

        line = 0
        for record in self.gb_records:
            for feature in record.features:
                if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                    if isinstance(feature.location, CompoundLocation):
                        if feature.strand == -1:
                            gene_start = int(feature.location.parts[-1].start.real)
                            gene_end = int(feature.location.parts[0].end.real)
                        else:
                            gene_start = int(feature.location.parts[0].start.real)
                            gene_end = int(feature.location.parts[-1].end.real)
                    else:
                        gene_start = int(feature.location.start.real)
                        gene_end = int(feature.location.end.real)

                    row_data = {
                        'Accession_Number': record.id,
                        'Sequence_length': len(record.seq),
                        'gene_start': gene_start,
                        'gene_end': gene_end,
                        'gene_strand': str(feature.strand),
                        'translation': feature.qualifiers['translation'][0],
                        'gene_order': line + 1,
                        'gene_id': feature.qualifiers.get('locus_tag', [''])[0],
                        'gene': feature.qualifiers.get('gene', ['-'])[0].lower() if 'gene' in feature.qualifiers else '-',
                        'product': feature.qualifiers.get('product', [''])[0]
                    }

                    for key, value in row_data.items():
                        genes_df.at[line, key] = value
                    line += 1

        return genes_df




    def solve_hypothetical_protein(self, genes_df):
        # 过滤假设蛋白
        genes_df['protein_length'] = genes_df['translation'].str.len()

        # 删除短的假设蛋白和短的非假设蛋白
        genes_df = genes_df[~((genes_df['product'] == 'hypothetical protein') & (genes_df['protein_length'] <= 50))]
        genes_df.reset_index(drop=True, inplace=True)
        # genes_df = genes_df[~((genes_df['product'] != 'hypothetical protein') & (genes_df['protein_length'] < 30))]
        # genes_df.reset_index(drop=True, inplace=True)

#############################################################################################
        # 新增：处理重叠的假设性蛋白，仅保留最长的（除非两者都>200）
        print("Handling overlapping hypothetical proteins...")
        
        # 获取所有假设性蛋白的索引
        hypo_indices = genes_df[genes_df['product'] == 'hypothetical protein'].index.tolist()
        
        # 预先标记哪些假设性蛋白需要删除（重叠情况）
        to_delete = set()
        
        # 检查每对假设性蛋白是否重叠
        for i, idx1 in enumerate(hypo_indices):
            if idx1 in to_delete:
                continue
                
            for idx2 in hypo_indices[i+1:]:
                if idx2 in to_delete:
                    continue
                    
                # 获取位置信息
                start1 = genes_df.at[idx1, 'gene_start']
                end1 = genes_df.at[idx1, 'gene_end']
                start2 = genes_df.at[idx2, 'gene_start']
                end2 = genes_df.at[idx2, 'gene_end']
                
                # 检查是否重叠
                if (start1 <= end2 and start2 <= end1):
                    # 获取蛋白长度
                    length1 = genes_df.at[idx1, 'protein_length']
                    length2 = genes_df.at[idx2, 'protein_length']
                    
                    # 如果两个假设性蛋白都大于200，则都保留
                    if length1 > 200 and length2 > 200:
                        print(f"Both hypothetical proteins at indices {idx1} and {idx2} overlap but both are >200aa, keeping both")
                        continue
                    
                    # 否则保留较长的蛋白
                    if length1 >= length2:
                        to_delete.add(idx2)
                        print(f"Marked hypothetical protein at index {idx2} for deletion (overlaps with {idx1}, shorter protein)")
                    else:
                        to_delete.add(idx1)
                        print(f"Marked hypothetical protein at index {idx1} for deletion (overlaps with {idx2}, shorter protein)")
                        # 如果第一个被标记删除，就不需要继续检查它与其他蛋白的重叠
                        break
        
        # 删除标记的重叠假设性蛋白
        if to_delete:
            genes_df = genes_df.drop(list(to_delete))
            genes_df.reset_index(drop=True, inplace=True)
            print(f"Removed {len(to_delete)} overlapping hypothetical proteins, keeping the longest ones")
################################################################################################################




        # 删除有overlap的假设蛋白
        genes_df['delete'] = 'no'
        index_list = genes_df[genes_df['product'] == 'hypothetical protein'].index

        # 添加调试信息
        print(f"Total genes after initial filtering: {len(genes_df)}")
        print(f"Indices of hypothetical proteins: {list(index_list)}")

        for i in index_list:
            # 检查 i + 1 是否在索引范围内
            if (i + 1) >= len(genes_df):
                print(f"Skipping index {i} as {i + 1} is out of range.")
                continue

            if (i != 0) and (i != genes_df.index[-1]):
                try:
                    end0 = genes_df.at[i-1, 'gene_end']
                    start1 = genes_df.at[i, 'gene_start']
                    end1 = genes_df.at[i, 'gene_end']
                    start2 = genes_df.at[i+1, 'gene_start']
                    if (start1 < end0) or (start2 < end1):
                        genes_df.at[i, 'delete'] = 'yes'
                except KeyError as e:
                    print(f"KeyError encountered at index {i}: {e}")
            elif i == 0:
                try:
                    start0 = genes_df.at[genes_df.index[-1], 'gene_start']
                    end0 = genes_df.at[genes_df.index[-1], 'gene_end']
                    start1 = genes_df.at[i, 'gene_start']
                    end1 = genes_df.at[i, 'gene_end']
                    start2 = genes_df.at[i + 1, 'gene_start']
                    if start0 > end0:
                        if (start1 < end0) or (start2 < end1):
                            genes_df.at[i, 'delete'] = 'yes'
                    elif start1 > end1:
                        if (start1 < end0) or (start2 < end1):
                            genes_df.at[i, 'delete'] = 'yes'
                    else:
                        if start2 < end1:
                            genes_df.at[i, 'delete'] = 'yes'
                except KeyError as e:
                    print(f"KeyError encountered at index {i} (first gene): {e}")
            elif i == genes_df.index[-1]:
                try:
                    end0 = genes_df.at[i-1, 'gene_end']
                    start1 = genes_df.at[i, 'gene_start']
                    end1 = genes_df.at[i, 'gene_end']
                    start2 = genes_df.at[0, 'gene_start']
                    end2 = genes_df.at[0, 'gene_end']
                    if start1 > end1:
                        if (start1 < end0) or (start2 < end1):
                            genes_df.at[i, 'delete'] = 'yes'
                    elif start2 > end2:
                        if (start1 < end0) or (start2 < end1):
                            genes_df.at[i, 'delete'] = 'yes'
                    else:
                        if start1 < end0:
                            genes_df.at[i, 'delete'] = 'yes'
                except KeyError as e:
                    print(f"KeyError encountered at index {i} (last gene): {e}")

        genes_df = genes_df.drop(genes_df[genes_df['delete'] == 'yes'].index)
        genes_df.reset_index(drop=True, inplace=True)

        # 删除RIP周围一定范围内的假设蛋白
        RIP_index_list = genes_df[genes_df['product'] == 'RIP'].index
        sequence_length = genes_df['Sequence_length'].iloc[0] if not genes_df.empty else 0
        indexes_to_drop = []

        for i in RIP_index_list:
            prev_indexes = [i - 2, i - 1]
            next_indexes = [i + 1, i + 2]

            for prev_index in prev_indexes:
                if prev_index >= 0:
                    try:
                        prev_product = genes_df.at[prev_index, 'product']
                        prev_distance = genes_df.at[i, 'gene_start'] - genes_df.at[prev_index, 'gene_end']
################################################################################
                        prev_gene_length = genes_df.at[prev_index,'protein_length']
################################################################################
                        if prev_product == 'hypothetical protein' and prev_distance < 100 and prev_gene_length < 150:
                            indexes_to_drop.append(prev_index)
                    except KeyError as e:
                        print(f"KeyError encountered when accessing previous index {prev_index}: {e}")
                else:
                    adjusted_prev_index = len(genes_df) + prev_index
                    if adjusted_prev_index >= 0:
                        try:
                            prev_product = genes_df.at[adjusted_prev_index, 'product']
                            prev_distance = (genes_df.at[i, 'gene_start'] + sequence_length - genes_df.at[adjusted_prev_index, 'gene_end'])
#################################################################################
                            prev_gene_length = genes_df.at[adjusted_prev_index, 'protein_length']
#################################################################################
                            if prev_product == 'hypothetical protein' and prev_distance < 100 and prev_gene_length < 150:
                                indexes_to_drop.append(adjusted_prev_index)
                        except KeyError as e:
                            print(f"KeyError encountered when accessing adjusted previous index {adjusted_prev_index}: {e}")

            for next_index in next_indexes:
                if next_index < len(genes_df):
                    try:
                        next_product = genes_df.at[next_index, 'product']
                        next_distance = genes_df.at[next_index, 'gene_start'] - genes_df.at[i, 'gene_end']
################################################################################
                        next_gene_length = genes_df.at[next_index,'protein_length']
################################################################################
                        if next_product == 'hypothetical protein' and next_distance < 100 and next_gene_length < 150:
                            indexes_to_drop.append(next_index)
                    except KeyError as e:
                        print(f"KeyError encountered when accessing next index {next_index}: {e}")
                else:
                    adjusted_next_index = next_index - len(genes_df)
                    if adjusted_next_index >= 0:
                        try:
                            next_product = genes_df.at[adjusted_next_index, 'product']
                            next_distance = (genes_df.at[adjusted_next_index, 'gene_start'] + sequence_length - genes_df.at[i, 'gene_end'])
#################################################################################
                            next_gene_length = genes_df.at[adjusted_next_index, 'protein_length']
#################################################################################
                            if next_product == 'hypothetical protein' and next_distance < 100 and next_gene_length < 150:
                                indexes_to_drop.append(adjusted_next_index)
                        except KeyError as e:
                            print(f"KeyError encountered when accessing adjusted next index {adjusted_next_index}: {e}")

        genes_df = genes_df.drop(indexes_to_drop)
        genes_df.reset_index(inplace=True, drop=True)
        solved_df = genes_df
        return solved_df

    def extract_intergenic_sequences(self, solved_df):
        intergenic_df = pd.DataFrame(columns=['Accession_Number', 'Intergenic_Start', 'Intergenic_End', 'Intergenic_Sequence'])
        gb_record = self.gb_records[0]
        seq_length = len(gb_record.seq)

        # 检查是否有足够的基因来提取 intergenic 序列
        if len(solved_df) < 2:
            # 如果没有基因或者只有一个基因，则将整个基因组作为一个 intergenic 序列
            intergenic_seq = str(gb_record.seq)
            if len(intergenic_seq) >= 50:
                new_row = {
                    'Accession_Number': solved_df.iloc[0]['Accession_Number'] if not solved_df.empty else gb_record.id,
                    'Intergenic_Start': 0,
                    'Intergenic_End': seq_length,
                    'Intergenic_Sequence': intergenic_seq
                }
                intergenic_df = pd.concat([intergenic_df, pd.DataFrame([new_row])], ignore_index=True)
            return intergenic_df

        # 如果有两个或更多基因，正常提取 intergenic 序列
        for i in range(len(solved_df) - 1):
            current_end = solved_df.iloc[i]['gene_end']
            next_start = solved_df.iloc[i+1]['gene_start']
            if current_end < next_start:
                intergenic_seq = str(gb_record.seq[current_end:next_start])
                if len(intergenic_seq) >= 50:
                    new_row = {
                        'Accession_Number': solved_df.iloc[i]['Accession_Number'],
                        'Intergenic_Start': current_end,
                        'Intergenic_End': next_start,
                        'Intergenic_Sequence': intergenic_seq
                    }
                    intergenic_df = pd.concat([intergenic_df, pd.DataFrame([new_row])], ignore_index=True)

        # 最后一个与第一个之间
        # last_end = solved_df.iloc[-1]['gene_end']
        # first_start = solved_df.iloc[0]['gene_start']
        # if last_end < first_start:
        #     intergenic_seq = str(gb_record.seq[last_end:first_start])
        # elif first_start < last_end and last_end > seq_length * 0.8:
        #     end_seq = str(gb_record.seq[last_end:])
        #     start_seq = str(gb_record.seq[:first_start])
        #     intergenic_seq = end_seq + start_seq
        # else:
        #     intergenic_seq = ""





########################################################################
        last_end = solved_df.iloc[-1]['gene_end']
        first_start = solved_df.iloc[0]['gene_start']
        seq_length = len(gb_record.seq)

        # 环状基因组中的基因间区判断
        if last_end < first_start:
            # 不跨越原点的基因间区
            intergenic_seq = str(gb_record.seq[last_end:first_start])
        else:
            # 两种可能：1.基因重叠 2.基因间区跨越原点
            # 计算"直接距离"和"环绕距离"
            direct_overlap = last_end - first_start  # 直接重叠的长度
            circular_distance = seq_length - last_end + first_start  # 环绕一周的距离
            
            if circular_distance > 0:  # 有环绕距离，表示有基因间区
                # 跨越原点的基因间区
                end_seq = str(gb_record.seq[last_end:])
                start_seq = str(gb_record.seq[:first_start])
                intergenic_seq = end_seq + start_seq
            else:
                # 真正的基因重叠，没有基因间序列
                intergenic_seq = ""
        if len(intergenic_seq) > 0.8 * seq_length:
            intergenic_seq = "" 
#########################################################################










        if len(intergenic_seq) >= 50:
            new_row = {
                'Accession_Number': solved_df.iloc[-1]['Accession_Number'],
                'Intergenic_Start': last_end,
                'Intergenic_End': first_start,
                'Intergenic_Sequence': intergenic_seq
            }
            intergenic_df = pd.concat([intergenic_df, pd.DataFrame([new_row])], ignore_index=True)

        # 加入RIP位置
        RIP_indices = solved_df[solved_df['product'] == 'RIP'].index
        for idx in RIP_indices:
            gene_start = solved_df.at[idx, 'gene_start']
            gene_end = solved_df.at[idx, 'gene_end']
            if gene_end < gene_start:
                sequence = str(gb_record.seq[gene_start:] + gb_record.seq[:gene_end])
            else:
                sequence = str(gb_record.seq[gene_start:gene_end])

            new_row = {
                'Accession_Number': solved_df.at[idx, 'Accession_Number'],
                'Intergenic_Start': gene_start,
                'Intergenic_End': gene_end,
                'Intergenic_Sequence': sequence
            }
            intergenic_df = pd.concat([intergenic_df, pd.DataFrame([new_row])], ignore_index=True)

        return intergenic_df

    def analyze_intergenic_sequences(self, intergenic_df):
        from_iteron = Iteron(intergenic_df, min_iteron=self.min_iteron, max_iteron=self.max_iteron)
        zcurve_analyzer = ZCurveAnalyzer(intergenic_df, sigma=self.sigma)
        pattern_finder = PatternFinder(intergenic_df)

        results = []
        for _, row in intergenic_df.iterrows():
            sequence = row['Intergenic_Sequence']

            iteron_score, iteron_sub_df, iteron_entropy_score = from_iteron.analyze(sequence)
            iteron_df_json = iteron_sub_df.to_json(orient='records')

            max_at_score, at_df = zcurve_analyzer.analyze(sequence)
            at_df_json = at_df.to_json(orient='records')

            pattern_score, RNA_df, pattern_sub_df = pattern_finder.find_patterns_in_sequence(sequence)
            RNA_df_json = RNA_df.to_json(orient='records')
            pattern_df_json = pattern_sub_df.to_json(orient='records')

            result_dict = {
                'Accession_Number': row['Accession_Number'],
                'Intergenic_Start': row['Intergenic_Start'],
                'Intergenic_End': row['Intergenic_End'],
                'Iteron_Score': iteron_score,
                'Iteron_Entropy_Score': iteron_entropy_score,
                'ZCurve_Score': max_at_score,
                'Pattern_Score': pattern_score,
                'iteron_df': iteron_df_json,
                'at_df': at_df_json,
                'RNA_df': RNA_df_json,
                'pattern_df': pattern_df_json
            }

            results.append(result_dict)

        merged_results_df = pd.DataFrame(results)
        return merged_results_df

    def merge_results(self, intergenic_df, merged_results_df):
        merged_df = intergenic_df.copy()
        merged_df = merged_df.merge(
            merged_results_df,
            on=['Accession_Number', 'Intergenic_Start', 'Intergenic_End'],
            how='left'
        )
        merged_df['Total_Score'] = (
            merged_df['Iteron_Score'].fillna(0) +
            merged_df['ZCurve_Score'].fillna(0) +
            abs(merged_df['Pattern_Score'].fillna(0))
        )
        return merged_df

    def run(self):
        genes_df = self.gbk2df()
        solved_df = self.solve_hypothetical_protein(genes_df)

        # 检查 solved_df 是否为空或行数不足
        if len(solved_df) < 1:
            print("erro")
        elif len(solved_df) < 2:
            print(f"erro {len(solved_df)} ")

        intergenic_df = self.extract_intergenic_sequences(solved_df)

        # 检查 intergenic_df 是否为空
        if intergenic_df.empty:

            gb_record = self.gb_records[0]
            seq_length = len(gb_record.seq)
            entire_seq = str(gb_record.seq)

            # 创建一个新的 intergenic 序列行
            new_row = {
                'Accession_Number': gb_record.id,
                'Intergenic_Start': 0,
                'Intergenic_End': seq_length,
                'Intergenic_Sequence': entire_seq
            }

            # 将新行添加到 intergenic_df
            intergenic_df = pd.concat([intergenic_df, pd.DataFrame([new_row])], ignore_index=True)

            print(intergenic_df.head())

        # 继续处理 intergenic_df 中的每个序列
        for index, row in intergenic_df.iterrows():
            gb_record = self.gb_records[0]
            intergenic_start = row['Intergenic_Start']
            intergenic_end = row['Intergenic_End']
            if intergenic_end < intergenic_start and len(solved_df) >= 2:
                intergenic_seq = gb_record.seq[intergenic_start:] + gb_record.seq[:intergenic_end]
            else:
                intergenic_seq = gb_record.seq[intergenic_start:intergenic_end]
            intergenic_df.at[index, 'Intergenic_Sequence'] = str(intergenic_seq)

        # **新增：删除长度大于等于3500的 intergenic_sequence**
        intergenic_df = intergenic_df[intergenic_df['Intergenic_Sequence'].apply(lambda s: 50 < len(s)) ]

        # 分析 intergenic 序列
        merged_results_df = self.analyze_intergenic_sequences(intergenic_df)

        # 合并结果
        merged_df = self.merge_results(intergenic_df, merged_results_df)

        return solved_df, merged_df



class OriginFinder:
    """
    主流程类：调用RIPPredictor和AnnotationAnalyzer对输入的FASTA文件进行分析，最终输出结果。
    """
    def __init__(self, fasta_file_path, db_file_path, hmm_file_path, output_dir,
                 e_value, min_score, min_alnlen,
                 min_iteron, max_iteron, sigma):
        self.fasta_file_path = fasta_file_path
        self.db_file_path = db_file_path
        self.hmm_file_path = hmm_file_path
        self.output_dir = output_dir
        

        # 运行RIPPredictor with optional parameters
        rip_predictor = RIPPredictor(
            fasta_file=self.fasta_file_path,
            db_file=self.db_file_path,
            hmm_file=self.hmm_file_path,
            output_dir=self.output_dir,
            e_value=e_value,
            min_score=min_score,
            min_alnlen=min_alnlen
        )
        self.gb_file_path, self.results_df_initial = rip_predictor.predict()

        self.sequence_length = self.get_sequence_length()
        self.min_iteron = min_iteron
        self.max_iteron = max_iteron
        self.sigma = sigma  # For ZCurveAnalyzer
        
        # 新增 OirV 数据库路径
        self.oriv_db_path = "OriVDB/OriV.fasta"

    def get_sequence_length(self):
        for record in SeqIO.parse(self.gb_file_path, 'genbank'):
            return len(record.seq)

    def get_gene_order_by_product(self, solved_df, product_name):
        gene_orders = solved_df[
            solved_df['product'].str.fullmatch(product_name, case=False, na=False)
        ]['gene_order'].tolist()
        return gene_orders if gene_orders else []

    def get_nearby_gene_orders(self, gene_order, range_value, max_gene_order):
        orders = []
        for i in range(-range_value, range_value + 1):
            nearby_order = gene_order + i
            if nearby_order < 1:
                nearby_order = max_gene_order + nearby_order
            elif nearby_order > max_gene_order:
                nearby_order = nearby_order - max_gene_order
            orders.append(nearby_order)
        return orders

    def calculate_distance(self, pos1_start, pos1_end, pos2_start, pos2_end):
        def normalize_position(pos):
            return pos % self.sequence_length

        pos1_start = normalize_position(pos1_start)
        pos1_end = normalize_position(pos1_end)
        pos2_start = normalize_position(pos2_start)
        pos2_end = normalize_position(pos2_end)

        def handle_crossing_origin(start, end):
            if end < start:
                return [(start, self.sequence_length), (0, end)]
            return [(start, end)]

        segments1 = handle_crossing_origin(pos1_start, pos1_end)
        segments2 = handle_crossing_origin(pos2_start, pos2_end)

        min_distance = float('inf')
        for seg1_start, seg1_end in segments1:
            for seg2_start, seg2_end in segments2:
                if seg1_end < seg2_start:
                    dist = seg2_start - seg1_end
                elif seg2_end < seg1_start:
                    dist = seg1_start - seg2_end
                else:
                    dist = 0
                min_distance = min(min_distance, dist)
        min_distance = min(min_distance,
                           self.sequence_length - max(pos1_end, pos2_end) + min(pos1_start, pos2_start))
        return min_distance
    
####################################################################增加blast
    def blast_intergenic_regions(self, merged_df):
        """
        对intergenic regions进行BLAST比对，保留最高匹配结果
        """
        import subprocess
        import tempfile
        import os
        from Bio import SeqIO
        import logging
        
        # 设置日志记录
        logging.basicConfig(
            filename=os.path.join(self.output_dir, "blast_debug.log"),
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        
        logging.info("开始BLAST比对分析，共有 %d 个intergenic regions", len(merged_df))
        
        # 从gbff文件中提取完整序列
        try:
            record = next(SeqIO.parse(self.gb_file_path, 'genbank'))
            genome_seq = record.seq
            logging.info("成功从GenBank文件读取序列，长度为: %d bp", len(genome_seq))
        except Exception as e:
            logging.error("读取GenBank文件出错: %s", str(e))
            return {}
        
        # 检查OirV数据库文件
        if not os.path.exists(self.oriv_db_path):
            logging.error("OirV数据库文件不存在: %s", self.oriv_db_path)
            return {}
        else:
            logging.info("OirV数据库文件位置: %s", self.oriv_db_path)
        
        # 创建临时文件用于BLAST
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_fasta:
            temp_fasta_path = temp_fasta.name
            logging.info("创建临时FASTA文件: %s", temp_fasta_path)
            
            # 为每个intergenic region创建序列
            seq_count = 0
            for idx, row in merged_df.iterrows():
                start = int(row['Intergenic_Start'])
                end = int(row['Intergenic_End'])
                
                # 处理跨越原点的情况
                if end < start:
                    seq = genome_seq[start:] + genome_seq[:end]
                    logging.debug("序列 %d 跨越原点: %d -> %d, 长度: %d bp", 
                                idx, start, end, len(seq))
                else:
                    seq = genome_seq[start:end]
                    logging.debug("序列 %d: %d -> %d, 长度: %d bp", 
                                idx, start, end, len(seq))
                
                temp_fasta.write(f">intergenic_{idx}\n{seq}\n")
                seq_count += 1
            
            logging.info("共写入 %d 个序列到临时FASTA文件", seq_count)
        
        # 创建临时文件用于BLAST输出
        temp_blast_out = tempfile.NamedTemporaryFile(delete=False).name
        logging.info("BLAST结果输出文件: %s", temp_blast_out)
        
        # 运行BLAST命令
        blast_cmd = [
            "blastn",
            "-query", temp_fasta_path,
            '-culling_limit', '1',
            "-subject", self.oriv_db_path,
            "-outfmt", "6 qseqid sseqid evalue pident qlen qstart qend slen length mismatch gapopen gaps bitscore",
            "-evalue", "1e-5",
            "-out", temp_blast_out
        ]
        
        logging.info("执行BLAST命令: %s", " ".join(blast_cmd))
        
        try:
            result = subprocess.run(blast_cmd, check=True, capture_output=True, text=True)
            if result.stderr:
                logging.warning("BLAST命令输出警告: %s", result.stderr)
            logging.info("BLAST命令执行成功")
        except subprocess.CalledProcessError as e:
            logging.error("BLAST命令执行失败: %s", str(e))
            logging.error("BLAST错误输出: %s", e.stderr)
            os.unlink(temp_fasta_path)
            if os.path.exists(temp_blast_out):
                os.unlink(temp_blast_out)
            return {}
        
        # 检查BLAST输出文件大小
        if os.path.getsize(temp_blast_out) == 0:
            logging.warning("BLAST输出文件为空，未找到任何匹配")
            os.unlink(temp_fasta_path)
            os.unlink(temp_blast_out)
            return {}
        
        # 解析BLAST结果
        blast_results = {}
        hit_count = 0
        filtered_count = 0
        
        with open(temp_blast_out, 'r') as blast_file:
            for line in blast_file:
                hit_count += 1
                cols = line.strip().split('\t')
                query_id = cols[0]
                idx = int(query_id.split('_')[1])
                
                coverage = float(cols[8]) / float(cols[7])
                identity = float(cols[3]) / 100.0
                evalue = float(cols[2])
                
                logging.debug("BLAST匹配 - 序列 %d: 匹配到 %s, E值: %g, 一致性: %.2f%%, 覆盖度: %.2f%%", 
                            idx, cols[1], evalue, identity*100, coverage*100)
                
                if coverage >= 0.75 and identity >= 0.8:
                    if idx not in blast_results or evalue < blast_results[idx]['evalue']:
                        blast_results[idx] = {
                            'hit_name': cols[1],
                            'evalue': evalue,
                            'identity': identity,
                            'coverage': coverage,
                            'start': int(cols[5]),
                            'end': int(cols[6])
                        }
                        logging.info("保留匹配 - 序列 %d: 匹配到 %s (E值: %g, 一致性: %.2f%%, 覆盖度: %.2f%%)", 
                                    idx, cols[1], evalue, identity*100, coverage*100)
                        filtered_count += 1
                else:
                    logging.debug("过滤掉不符合要求的匹配 (coverage < 0.75 或 identity < 0.8)")
        
        logging.info("BLAST匹配总数: %d, 满足条件的匹配数: %d", hit_count, filtered_count)
        
        # 清理临时文件
        os.unlink(temp_fasta_path)
        os.unlink(temp_blast_out)
        logging.info("清理临时文件完成")
        
        return blast_results



    def analyze_intergenic_regions(self, solved_df, merged_df):
        # 新增：对intergenic regions进行BLAST比对
        blast_results = self.blast_intergenic_regions(merged_df)

        # 初始化结果列表和RIP信息
        results = []
        temp_results = {}  # 临时存储评分结果
        max_gene_order = solved_df['gene_order'].max()

        # 提取RIP和分区蛋白相关信息
        rip_orders = self.get_gene_order_by_product(solved_df, 'RIP')
        partition_pattern = r'partition protein|ParA|ParB'
        partitioning_orders = solved_df[solved_df['product'].str.contains(partition_pattern, case=False, na=False, regex=True)]['gene_order'].tolist()

        # 提取RIP位置
        rip_positions = []
        for _, gene in solved_df.iterrows():
            if gene['gene_order'] in rip_orders:
                rip_positions.append((gene['gene_start'], gene['gene_end']))

        # 构建RIP相邻基因关联
        rip_nearby_orders = []
        rip_association = {}
        # 新增：构建RIP与nearby_RIP的关联映射
        nearby_regions_by_rip = {rip_order: [] for rip_order in rip_orders}  # 存储每个RIP关联的nearby区域
        
        for rip_order in rip_orders:
            nearby_orders = self.get_nearby_gene_orders(rip_order, 1, max_gene_order)
            for o in nearby_orders:
                if o not in rip_association:
                    rip_association[o] = rip_order
            rip_nearby_orders.extend(nearby_orders)

        # 其他初始化代码
        rip_nearby_orders = list(set(rip_nearby_orders))
        partitioning_nearby_orders = []
        for part_order in partitioning_orders:
            nearby_orders = self.get_nearby_gene_orders(part_order, 6, max_gene_order)
            partitioning_nearby_orders.extend(nearby_orders)
        partitioning_nearby_orders = list(set(partitioning_nearby_orders))
        
        # First pass for partitioning regions
        partitioning_regions = []
        for _, intergenic in merged_df.iterrows():
            start = intergenic['Intergenic_Start']
            end = intergenic['Intergenic_End']
            pattern_score = intergenic['Pattern_Score']

            for _, gene in solved_df.iterrows():
                gene_start = gene['gene_start']
                gene_end = gene['gene_end']
                gene_order = gene['gene_order']

                min_distance = self.calculate_distance(gene_start, gene_end, start, end)
                if min_distance >= 10:
                    continue
                if gene_order in partitioning_nearby_orders:
                    partitioning_regions.append({
                        'start': start,
                        'end': end,
                        'pattern_score': pattern_score,
                        'total_score': intergenic['Total_Score']
                    })
        
        max_score_with_pattern = None
        max_score_without_pattern = None
        if partitioning_regions:
            pattern_regions = [r for r in partitioning_regions if abs(r['pattern_score']) >= 10]
            if pattern_regions:
                max_score_with_pattern = max(pattern_regions, key=lambda x: x['total_score'])
            max_score_without_pattern = max(partitioning_regions, key=lambda x: x['total_score'])

        # ========== 第一次循环：处理所有非RIP区域 ==========
        rip_regions = {}  # 存储RIP区域信息，用于第二次循环
        
        for idx, intergenic in merged_df.iterrows():
            start = intergenic['Intergenic_Start']
            end = intergenic['Intergenic_End']
            level = 0
            evidence = "N/A"
            pattern_score = intergenic['Pattern_Score']
            iteron_score = intergenic['Iteron_Score']
            at_score = intergenic['ZCurve_Score']

            # 检查是否为RIP区域
            is_rip_region = any((start == rip_start and end == rip_end) for rip_start, rip_end in rip_positions)
            associated_rip_order = None
            
            # BLAST匹配检查
            blast_hit = None
            if idx in blast_results:
                blast_hit = blast_results[idx]
                level = max(level, 4)
                evidence = "blast_hit"
                
            if is_rip_region:
                # 先将RIP区域存储起来，稍后处理
                matching_gene = solved_df[(solved_df['gene_start'] == start) & (solved_df['gene_end'] == end)]
                if not matching_gene.empty:
                    associated_rip_order = matching_gene.iloc[0]['gene_order']
                    rip_regions[associated_rip_order] = {
                        'idx': idx,
                        'intergenic': intergenic,
                        'start': start,
                        'end': end,
                        'blast_hit': blast_hit
                    }
            else:
                # 非RIP区域的处理
                is_nearby_rip = False
                is_nearby_par = False

                for _, gene in solved_df.iterrows():
                    gene_start = gene['gene_start']
                    gene_end = gene['gene_end']
                    gene_order = gene['gene_order']

                    min_distance = self.calculate_distance(gene_start, gene_end, start, end)
                    if min_distance >= 10:
                        continue

                    if gene_order in rip_orders or gene_order in rip_nearby_orders:
                        if gene_order in rip_orders:
                            associated_rip_order = gene_order
                        elif gene_order in rip_nearby_orders and gene_order in rip_association:
                            associated_rip_order = rip_association[gene_order]

                        # 如果是nearby_RIP并且满足条件，设置level为4
                        if (iteron_score >= 20 or abs(pattern_score) >= 10 or
                                pattern_score < 0 or pattern_score >= 1000 or
                                at_score >= 20):
                            level = max(level, 4)
                            is_nearby_rip = True
                            
                            # 新增：记录与RIP关联的区域信息
                            if associated_rip_order in nearby_regions_by_rip:
                                nearby_regions_by_rip[associated_rip_order].append({
                                    'idx': idx,
                                    'level': level,
                                    'evidence': "nearby RIP" if not blast_hit else evidence
                                })
                                
                    elif gene_order in partitioning_nearby_orders:
                        # 分区蛋白相关逻辑保持不变
                        if pattern_score >= 1000:
                            level = max(level, 4)
                            is_nearby_par = True
                        elif (max_score_with_pattern and
                            start == max_score_with_pattern['start'] and
                            end == max_score_with_pattern['end']):
                            level = max(level, 3)
                            is_nearby_par = True
                        elif (max_score_without_pattern and
                            start == max_score_without_pattern['start'] and
                            end == max_score_without_pattern['end']):
                            level = max(level, 3)
                            is_nearby_par = True

                if is_nearby_rip and not blast_hit:
                    evidence = "nearby RIP"
                elif is_nearby_par and not blast_hit:
                    evidence = "nearby PAR"

                # 特殊情况：高模式得分
                if pattern_score >= 1000 and not blast_hit:
                    evidence = "RNAII"
                    level = max(level, 4)

                # 更新总分和其他信息
                if blast_hit:
                    total_score = 1000
                else:
                    if level == 0:
                        sum_score = (max(intergenic['Pattern_Score'], 0) +
                                    intergenic['Iteron_Score'] +
                                    intergenic['ZCurve_Score'])
                    else:
                        sum_score = 0
                    total_score = intergenic['Total_Score']

                # 创建结果字典
                result = intergenic.to_dict()
                result.update({
                    'Level': level,
                    'Evidence': evidence,
                    'Sum_Score': sum_score if 'sum_score' in locals() else 0,
                    'Associated_RIP_Order': associated_rip_order,
                    'Total_Score': total_score
                })
                
                # 添加BLAST匹配信息
                if blast_hit:
                    result.update({
                        'blast_hit_name': blast_hit['hit_name'],
                        'blast_evalue': blast_hit['evalue'],
                        'blast_identity': blast_hit['identity'],
                        'blast_coverage': blast_hit['coverage'],
                        'blast_start': blast_hit['start'],
                        'blast_end': blast_hit['end']
                    })
                
                # 存储非RIP区域的结果
                temp_results[idx] = result

        # ========== 第二次循环：处理RIP区域 ==========
        for rip_order, rip_info in rip_regions.items():
            idx = rip_info['idx']
            intergenic = rip_info['intergenic']
            start = rip_info['start']
            end = rip_info['end']
            blast_hit = rip_info['blast_hit']
            
            # 默认不将RIP设为Level 4
            level = 0  
            evidence = "RIP"
            
            # 检查该RIP关联的nearby_RIP区域是否有Level 4
            nearby_level_4_exists = False
            if rip_order in nearby_regions_by_rip:
                for nearby in nearby_regions_by_rip[rip_order]:
                    if nearby['level'] == 4:
                        nearby_level_4_exists = True
                        break
                        
            # 只有在没有Level 4的nearby_RIP时，才将RIP设为Level 4
            if not nearby_level_4_exists or blast_hit:
                level = max(level, 4)
                
            # 如果有BLAST匹配，优先使用BLAST证据
            if blast_hit:
                evidence = "blast_hit"
                total_score = 1000
            else:
                pattern_score = intergenic['Pattern_Score']
                iteron_score = intergenic['Iteron_Score']
                at_score = intergenic['ZCurve_Score']
                
                if level == 0:
                    sum_score = (max(pattern_score, 0) + iteron_score + at_score)
                else:
                    sum_score = 0
                total_score = intergenic['Total_Score']
            
            # 创建结果字典
            result = intergenic.to_dict()
            result.update({
                'Level': level,
                'Evidence': evidence,
                'Sum_Score': sum_score if 'sum_score' in locals() else 0,
                'Associated_RIP_Order': rip_order,
                'Total_Score': total_score
            })
            
            # 添加BLAST匹配信息
            if blast_hit:
                result.update({
                    'blast_hit_name': blast_hit['hit_name'],
                    'blast_evalue': blast_hit['evalue'],
                    'blast_identity': blast_hit['identity'],
                    'blast_coverage': blast_hit['coverage'],
                    'blast_start': blast_hit['start'],
                    'blast_end': blast_hit['end']
                })
            
            # 存储RIP区域的结果
            temp_results[idx] = result
        
        # 将所有结果按原始顺序添加到results列表
        for idx in sorted(temp_results.keys()):
            results.append(temp_results[idx])
        
        # 构建结果DataFrame并整理列顺序
        results_df = pd.DataFrame(results)
        columns = list(results_df.columns)
        intergenic_end_index = columns.index('Intergenic_End')
        if 'Level' in columns:
            columns.remove('Level')
        if 'Evidence' in columns:
            columns.remove('Evidence')
        columns.insert(intergenic_end_index + 1, 'Level')
        columns.insert(intergenic_end_index + 2, 'Evidence')
        
        # 重新排列BLAST相关列
        blast_columns = [col for col in columns if col.startswith('blast_')]
        for col in blast_columns:
            columns.remove(col)
        columns.extend(blast_columns)
        
        results_df = results_df[columns]

        # 分级过滤
        level_0_mask = results_df['Level'] == 0
        level_0_df = results_df[level_0_mask].copy()
        other_levels_df = results_df[~level_0_mask].copy()

        other_levels_df = other_levels_df.sort_values(['Level', 'Total_Score'], ascending=[False, False])
        other_levels_df['Sum_Score'] = 0

        level_4_df = other_levels_df[other_levels_df['Level'] == 4].copy()
        non_level_4_df = other_levels_df[other_levels_df['Level'] != 4]

        filter_evidences = ["RIP", "nearby RIP"]
        level_4_need_filter = level_4_df[level_4_df['Evidence'].isin(filter_evidences)]
        level_4_no_filter = level_4_df[~level_4_df['Evidence'].isin(filter_evidences)]

        # 修改：优先保留nearby_RIP区域
        def keep_max_score_per_rip_order(group):
            # 首先检查组内是否有 "nearby RIP" 证据的区域
            nearby_rip_rows = group[group['Evidence'] == "nearby RIP"]
            
            if not nearby_rip_rows.empty:
                # 如果有 "nearby RIP" 证据的区域，选择其中总分最高的
                max_nearby_score = nearby_rip_rows['Total_Score'].max()
                return nearby_rip_rows[nearby_rip_rows['Total_Score'] == max_nearby_score]
            else:
                # 如果没有 "nearby RIP" 证据的区域，则从整个组中选择总分最高的
                max_score = group['Total_Score'].max()
                return group[group['Total_Score'] == max_score]

        level_4_with_rip_order = level_4_need_filter[level_4_need_filter['Associated_RIP_Order'].notna()]
        level_4_without_rip_order = level_4_need_filter[level_4_need_filter['Associated_RIP_Order'].isna()]

        level_4_with_rip_order = level_4_with_rip_order.groupby('Associated_RIP_Order', group_keys=False).apply(keep_max_score_per_rip_order, include_groups=False)

        level_4_filtered = pd.concat([level_4_with_rip_order, level_4_without_rip_order, level_4_no_filter], ignore_index=True)
        level_4_filtered = level_4_filtered.sort_values(['Level', 'Total_Score'], ascending=[False, False])

        other_levels_df = pd.concat([non_level_4_df, level_4_filtered], ignore_index=True)
        other_levels_df = other_levels_df.sort_values(['Level', 'Total_Score'], ascending=[False, False])

        level_0_df['Pattern_Score'] = level_0_df['Pattern_Score'].clip(lower=0)
        level_0_df['Sum_Score'] = (level_0_df['Pattern_Score'] +
                                level_0_df['Iteron_Score'] +
                                level_0_df['ZCurve_Score'])
        level_0_df = level_0_df.sort_values('Sum_Score', ascending=False)

        results_df = pd.concat([other_levels_df, level_0_df])
        # results_df.to_csv(os.path.join(self.output_dir, "test.csv"), index=False)
        
        # 过滤低分长区域
        intergenic_length = results_df['Intergenic_End'] - results_df['Intergenic_Start']
        results_df = results_df[~((results_df['Total_Score'] < 500) & (intergenic_length > 1500))]
        
        # 检查Level分布并选择最终结果
        level_counts = results_df['Level'].value_counts()

        if 4 in level_counts and level_counts[4] > 0:
            # 存在 Level==4 时（对应 Type 1）只保留 Level==4 的行
            results_df = results_df[results_df['Level'] == 4].copy()
        elif 3 in level_counts and level_counts[3] > 0:
            # 不存在 Level==4 但存在 Level==3 时（对应 Type 2），对 Level==3 重新计算 Sum_Score，
            # 并保留 Sum_Score 最大的行
            level_3_df = results_df[results_df['Level'] == 3].copy()
            level_3_df['Sum_Score'] = (np.maximum(level_3_df['Pattern_Score'], 0) 
                                    + level_3_df['Iteron_Score'] 
                                    + level_3_df['ZCurve_Score'])
            max_sum_3 = level_3_df['Sum_Score'].max()
            results_df = level_3_df[level_3_df['Sum_Score'] == max_sum_3].copy()
        else:
            # 不存在 Level==4 和 Level==3 时，只保留 Level==0 中 Sum_Score 最大的行（对应 Type 3）
            level_0_df = results_df[results_df['Level'] == 0].copy()
            max_sum_0 = level_0_df['Sum_Score'].max()
            results_df = level_0_df[level_0_df['Sum_Score'] == max_sum_0].copy()

        # 将 Level 映射为 Type
        results_df['Type'] = results_df['Level'].map({4: "1", 3: "2", 0: "3"})
        results_df = results_df.drop(columns=['Level'])
#########################################################################################################################################
        # 创建包含所有结果的副本并分配Type
        full_results_df = pd.concat([other_levels_df, level_0_df])
        full_results_df['Type'] = full_results_df['Level'].map({4: "1", 3: "2", 0: "3"})

        # 删除Level列和Iteron_Score列
        full_results_df = full_results_df.drop(columns=['Level', 'Iteron_Score', 'Associated_RIP_Order'])

        # 重命名列
        rename_map = {
            'iteron_df': 'Iteron_df',
            'pattern_df': 'Pattern_Df',
            'ZCurve_Score': 'AT_Score',
            'at_df': 'AT_Df',
            'blast_hit_name': 'Blast_Hit',
            'blast_evalue': 'Blast_Evalue',
            'blast_identity': 'Blast_Identity',
            'blast_coverage': 'Blast_Coverage',
            'blast_start': 'Blast_Start',
            'blast_end': 'Blast_End'
        }
        full_results_df = full_results_df.rename(columns=rename_map)

        # 定义列的排序
        ordered_columns = [
            'Accession_Number', 'Evidence', 'Type', 'Intergenic_Start', 
            'Intergenic_End', 'Intergenic_Sequence', 'Iteron_Entropy_Score', 
            'Iteron_df', 'Pattern_Score', 'Pattern_Df', 'AT_Score', 'AT_Df',
            'Blast_Hit', 'Blast_Evalue', 'Blast_Identity', 'Blast_Coverage',
            'Blast_Start', 'Blast_End'
        ]

        # 确保所有列都存在，不存在的列忽略
        available_columns = [col for col in ordered_columns if col in full_results_df.columns]
        other_columns = [col for col in full_results_df.columns if col not in ordered_columns]
        final_columns = available_columns + other_columns

        # 重新排序列
        full_results_df = full_results_df[final_columns]

        # 按照Type从1到3排列，同一Type内按照Total_Score降序排列
        full_results_df = full_results_df.sort_values(by=['Type', 'Total_Score'], ascending=[True, False])

        # 将所有带Type的intergenic结果输出为CSV
        full_results_df.to_csv(os.path.join(self.output_dir, "All_IGSs.csv"), index=False)

        # 最终结果按Type和Total_Score排序
        results_df = results_df.sort_values(by=['Type', 'Total_Score'], ascending=[True, False])
        # results_df.to_csv(os.path.join(self.output_dir, "selected_ori_regions.csv"), index=False)
###########################################################################################################################################





        # 创建位置函数
        def create_location(start, end, seq_length, strand=+1):
            if start < 0 or end < 0:
                raise ValueError("Start and end must be non-negative.")
            
            # 确保位置在序列长度范围内
            start = start % seq_length
            end = end % seq_length
            
            # 如果end == start表明覆盖整条环
            if start == end:
                return FeatureLocation(0, seq_length, strand=strand)

            if end < start:
                # 如果end == 0，则直接返回从start到序列末尾的FeatureLocation
                if end == 0:
                    return FeatureLocation(start, seq_length, strand=strand)
                else:
                    # 正常跨越原点处理
                    # 注意：这里不要使用seq_length作为part1的终点，而是使用seq_length-1
                    part1 = FeatureLocation(start, seq_length-1, strand=strand)
                    part2 = FeatureLocation(0, end, strand=strand)
                    return CompoundLocation([part1, part2])
            else:
                return FeatureLocation(start, end, strand=strand)

        # 读取GenBank文件并注释
        records = list(SeqIO.parse(self.gb_file_path, "genbank"))
        if len(records) != 1:
            raise ValueError("Currently only handling single-record GBFF files.")
        record = records[0]
        seq_length = len(record.seq)

        # 对最终过滤后的results_df逐行处理进行注释
        for _, row in results_df.iterrows():
            intergenic_start = int(row['Intergenic_Start'])
            intergenic_end = int(row['Intergenic_End'])

            # 注释oriV为misc_feature，note中包含Type和Evidence
            oriV_feature = SeqFeature(
                location=create_location(intergenic_start, intergenic_end, seq_length, strand=1),
                type="misc_feature"
            )
            
            # 添加注释信息
            note_text = f"OriV; Type={row['Type']}; Evidence={row['Evidence']}"
            if 'blast_hit_name' in row and pd.notna(row['blast_hit_name']):
                note_text += f"; BLAST_hit={row['blast_hit_name']}; E-value={row['blast_evalue']}; Identity={row['blast_identity']:.2f}; Coverage={row['blast_coverage']:.2f}"
            
            oriV_feature.qualifiers['note'] = [note_text]
            record.features.append(oriV_feature)

            # 解析和注释其他特征
            iteron_data = json.loads(row['iteron_df']) if 'iteron_df' in row and pd.notna(row['iteron_df']) else []
            at_data = json.loads(row['at_df']) if 'at_df' in row and pd.notna(row['at_df']) else []
            RNA_data = json.loads(row['RNA_df']) if 'RNA_df' in row and pd.notna(row['RNA_df']) else []
            pattern_data = json.loads(row['pattern_df']) if 'pattern_df' in row and pd.notna(row['pattern_df']) else []

            # 注释iteron
            for it in iteron_data:
                rel_start = int(it['start'])    # 1-based起点
                rel_end = int(it['end'])        # 1-based终点
                # 转换为0-based, end为排除端：起点减1，终点不减1
                abs_start = (intergenic_start + rel_start - 1) % seq_length
                abs_end = (intergenic_start + rel_end) % seq_length
                it_feature = SeqFeature(
                    location=create_location(abs_start, abs_end, seq_length, strand=1),
                    type="misc_feature"
                )
                it_feature.qualifiers['note'] = ["hypothetical_iteron"]
                record.features.append(it_feature)

            # 注释AT富集区
            for at in at_data:
                rel_start = int(at['start_n'])
                rel_end = int(at['end_n'])
                abs_start = (intergenic_start + rel_start - 1) % seq_length
                abs_end = (intergenic_start + rel_end) % seq_length
                at_feature = SeqFeature(
                    location=create_location(abs_start, abs_end, seq_length, strand=1),
                    type="misc_feature"
                )
                at_feature.qualifiers['note'] = ["AT_rich_region"]
                record.features.append(at_feature)

            # 注释RNA结构
            for rna in RNA_data:
                rel_start = int(rna['start'])
                rel_end = int(rna['end'])
                abs_start = (intergenic_start + rel_start - 1) % seq_length
                abs_end = (intergenic_start + rel_end) % seq_length
                rna_strand = 1 if rna['strand'] == '+' else -1
                rna_feature = SeqFeature(
                    location=create_location(abs_start, abs_end, seq_length, strand=rna_strand),
                    type="misc_feature"
                )
                # 分别存放名称和E-value
                rna_feature.qualifiers['note'] = [rna['name']]
                rna_feature.qualifiers['e-value'] = [rna['evalue']]
                record.features.append(rna_feature)

            # 注释模式
            for pat in pattern_data:
                rel_start = int(pat['start'])
                rel_end = int(pat['end'])
                abs_start = (intergenic_start + rel_start - 1) % seq_length
                abs_end = (intergenic_start + rel_end) % seq_length
                pat_strand = 1 if pat['strand'] == '+' else -1
                pat_feature = SeqFeature(
                    location=create_location(abs_start, abs_end, seq_length, strand=pat_strand),
                    type="misc_feature"
                )
                # 仅存放名称，无需e-value
                pat_feature.qualifiers['note'] = [pat['name']]
                record.features.append(pat_feature)
                
            # 注释BLAST匹配
            if 'blast_hit_name' in row and pd.notna(row['blast_hit_name']):
                rel_start = int(row['blast_start']) if pd.notna(row['blast_start']) else 1
                rel_end = int(row['blast_end']) if pd.notna(row['blast_end']) else rel_start + 100  # 默认长度
                abs_start = (intergenic_start + rel_start - 1) % seq_length
                abs_end = (intergenic_start + rel_end) % seq_length
                
                blast_feature = SeqFeature(
                    location=create_location(abs_start, abs_end, seq_length, strand=1),
                    type="misc_feature"
                )
                blast_feature.qualifiers['note'] = [f"BLAST_hit: {row['blast_hit_name']}"]
                blast_feature.qualifiers['e-value'] = [str(row['blast_evalue'])]
                blast_feature.qualifiers['identity'] = [f"{row['blast_identity']:.2f}"]
                blast_feature.qualifiers['coverage'] = [f"{row['blast_coverage']:.2f}"]
                record.features.append(blast_feature)

        # 将修改后的记录写回gbff文件
        with open(self.gb_file_path, 'w') as output_handle:
            SeqIO.write(record, output_handle, 'genbank')

        return results_df





















    def find_origin(self):
        annotation_analyzer = AnnotationAnalyzer(self.gb_file_path,min_iteron=self.min_iteron,max_iteron=self.max_iteron,sigma=self.sigma)
        solved_df, merged_df = annotation_analyzer.run()
        
        results_df = self.analyze_intergenic_regions(solved_df, merged_df)


        # 打印统计信息
        print("\nType Distribution:")
        print(results_df['Type'].value_counts().sort_index())

        print("\nHighest Type Sequences:")
        highest_type = results_df['Type'].max()
        print(results_df[results_df['Type'] == highest_type][
            ['Intergenic_Start', 'Intergenic_End', 'Type', 'Evidence', 'Total_Score']
        ].to_string())

        # 将更新后的gbff文件复制到results目录下
        # self.gb_file_path是更新过RIP标记的gbff文件
        gbff_copy_path = os.path.join(self.output_dir, os.path.basename(self.gb_file_path))
        shutil.copyfile(self.gb_file_path, gbff_copy_path)
        print(f"Updated gbff file copied to {gbff_copy_path}")

        return solved_df, merged_df, results_df










############################################################
#                        Main
############################################################
    

    
#oriVfinder.py
def main():
    parser = argparse.ArgumentParser(description="Run OriginFinder pipeline")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("--db_file", default="../orivfinder/app/RIPDB/RIPDB", help="Path to MMseqs2 database file")
    parser.add_argument("--hmm_file", default="../orivfinder/app/hmmdb/912.hmm", help="Path to HMM profile file")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory")


    # Additional optional arguments for OriginFinder
    parser.add_argument("--e_value", type=float, default=1e-5, help="E-value threshold for RIPPredictor")
    parser.add_argument("--min_score", type=float, default=50, help="Minimum score for RIPPredictor")
    parser.add_argument("--min_alnlen", type=int, default=100, help="Minimum alignment length for RIPPredictor")
    parser.add_argument("--min_iteron", type=int, default=11, help="Minimum iteron length for Iteron")
    parser.add_argument("--max_iteron", type=int, default=22, help="Maximum iteron length for Iteron")
    parser.add_argument("--sigma", type=float, default=5, help="Sigma value for ZCurveAnalyzer")

    args = parser.parse_args()

    origin_finder = OriginFinder(
        fasta_file_path=args.fasta,
        db_file_path=args.db_file,
        hmm_file_path=args.hmm_file,
        output_dir=args.output_dir,
        e_value=args.e_value,
        min_score=args.min_score,
        min_alnlen=args.min_alnlen,
        min_iteron=args.min_iteron,
        max_iteron=args.max_iteron,
        sigma=args.sigma
    )

    solved_df, merged_df, results_df = origin_finder.find_origin()

    print("Processing completed. Results are saved in:", args.output_dir)


if __name__ == "__main__":
    main()
