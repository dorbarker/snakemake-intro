import subprocess
import sys
import pandas as pd
import os
from pathlib import Path
import itertools
# Add the scripts subdirectory to $PATH
scripts = str(Path(workflow.current_basedir).joinpath("scripts").resolve()) + os.pathsep
os.environ["PATH"] = scripts + os.environ["PATH"]

def get_version(shell_command: str) -> str:
	ret = subprocess.check_output(shell_command, shell = True)
	s = ret.decode().replace("\n", " - ") # no newlines in version strings
	return s
# Versions
PYTHON_VERSION = sys.version.split()[0]
PANDAS_VERSION = pd.__version__


# Shared data
seeds = list(range(config["balanced_cohort_number"]))
primer_sizes = list(range(150, 1000, 100))

rule all:
	input:
		"results/ntb.tsv",
		"results/positive_hits.tsv",
		"results/primer_pair_calls.tsv"

rule format_grapetree_calls:
	input:
		"calls.tsv"
	output:
		"grapetree_calls.tsv"
	threads: 1
	params:
		mem="1G",
		time="01:00"
	shell:
		'vim -E	+"norm! ggi#" +"w {output}" -scq! {input}'

rule calculate_newick_tree:
	input:
		rules.format_grapetree_calls.output
	output:
		"tree.nwk"
	threads:
		5
	params:
		mem="16G",
		time="02:00:00"
	conda: "envs/grapetree.yaml"
	shell:
		"grapetree --profile {input} > {output}"

rule cluster_grapetree:
	input:
		rules.calculate_newick_tree.output
	output:
		"clusters.tsv"
	threads:
		1
	params:
		mem="4G",
		time="06:00:00"
	shell:
		"gtclust -i {input} -o {output}"

rule calculate_cohorts:
	input:
		clusters=rules.cluster_grapetree.output,
		metadata="metadata.tsv"
	output:
		"tabular_cohorts.tsv"
	threads:
		1
	params:
		mem="16G",
		time="04:00:00",
		size_1=config["size_1"],
		size_2=config["size_2"],
		pos_1=config["pos_1"],
		pos_2=config["pos_2"],
		neg_1=config["neg_1"],
		neg_2=config["neg_2"],
		variable=config["variable"],
		posval=config["positive_value"]
	shell:
		"cohortid "
		"--clusters {input.clusters} "
		"--metadata {input.metadata} "
		"--output {output} "
		"--variable {params.variable} "
		"--positive-value {params.posval} "
		"--size-1 {params.size_1} --size-2 {params.size_2} "
		"--pos-1 {params.pos_1} --pos-2 {params.pos_2} "
		"--neg-1 {params.neg_1} --neg-2 {params.neg_2} "

rule expand_cohorts:
	input:
		clusters=rules.cluster_grapetree.output,
		cohorts=rules.calculate_cohorts.output
	output:
		"cohorts.tsv"
	params:
		mem="4G",
		time="10:00"
	conda:
		"envs/r.yaml"
	shell:
		"cohort-expander.R "
		"--clusters {input.clusters} "
		"--cohorts {input.cohorts} "
		"--output {output}"

rule pirate_pangenome:
	output:
		"pangenome/PIRATE.gene_families.ordered.tsv"

	threads:
		72
	params:
		time="5-0",
		mem="144G"
	conda: "envs/pirate.yaml"
	shell:
		"PIRATE "
		"--input gffs/ "
		"--output pangenome/ "
		"--nucl "
		"--align "
		"--para-off "
		"--threads {threads} "

rule get_PIRATE_to_Rtab:
	output: "scripts/PIRATE_to_Rtab.pl"
	threads: 1
	params:
		mem="100M",
		time="10:00",
		url="SionBayliss/PIRATE/master/tools/convert_format/PIRATE_to_Rtab.pl"
	shell: "curl -o {output} https://raw.githubusercontent.com/{params.url}"

rule binarize_markers:
	input:
		script=rules.get_PIRATE_to_Rtab.output,
		pirate_results="pangenome/PIRATE.gene_families.ordered.tsv"
	output:
		"feht_inputs/binary_calls.tsv"
	threads: 1
	params:
		time="30:00",
		mem="8G"
	shell:
		"{input.script} -i pangenome/PIRATE.*.tsv -o {output}"

rule transpose_binary_markers:
	input:
		rules.binarize_markers.output
	output:
		"genes_binarized.tsv"
	threads: 1
	params:
		mem="4G",
		time="01:00:00"
	run:
		import pandas as pd
		original = pd.read_csv(input[0], sep="\t", header=0)
		transposed = original.transpose()
		transposed.to_csv(output[0], sep="\t", header=False)

rule find_maximum_negative_threshold:
	input:
		rules.expand_cohorts.output
	output:
		"feht_inputs/maxthreshold.txt"
	params:
		mem="4G",
		time="10:00"
	run:
		import pandas as pd
		from pathlib import Path
		cohorts = pd.read_csv(input[0], sep = "\t", header=0, index_col=0)
		#pos = cohorts.apply(lambda col: col.eq("positive").sum())
		neg = cohorts.apply(lambda col: col.eq("negative").sum())
		#most_balanced = abs(((neg / pos) - 1)).values.argmin()
		#most_positive_colname = pos.index[pos.values.argmax()]
		most_negative_colname = neg.index[neg.values.argmax()]
		#most_balanced_colname = neg.index[most_balanced]
		#Path(output[0]).write_text(most_balanced_colname)
		#Path(output[0]).write_text(most_positive_colname)
		Path(output[0]).write_text(most_negative_colname)

rule run_feht:
	input:
		maxthresh=rules.find_maximum_negative_threshold.output,
		markers=rules.binarize_markers.output,
		cohorts=rules.expand_cohorts.output
	output:
		"feht_results.tsv"
	threads: 1
	params:
		#threshold=storage.fetch("max_negative_threshold"),
		mem="16G",
		time="2-0"
	conda: "envs/feht.yaml"
	shell:
		"mapfile -t < {input.maxthresh} maxthresh && "
		"feht "
		"-i {input.cohorts} "
		"-d {input.markers} "
		"--one \"${{maxthresh[0]}} positive none\" "
		"--two \"${{maxthresh[0]}} negative\" "
		"--mode binary "
		"--correction bonferroni "
		"> {output} "

rule calculate_ntb:
	input:
		rules.run_feht.output
	output:
		"results/ntb.tsv"
	params:
		mem="4G",
		time="30:00"
	shell:
		"ntb-feht.py --input {input} --output {output}"

rule create_balanced_cohorts:
	input:
		cohorts=rules.expand_cohorts.output,
		maxthresh=rules.find_maximum_negative_threshold.output,
		calls=rules.transpose_binary_markers.output,
		metadata="metadata.tsv"
	output:
		metadata="balanced_cohorts/{seed}/metadata.tsv",
		calls="balanced_cohorts/{seed}/calls.tsv",
		cohorts="balanced_cohorts/{seed}/cohorts.tsv"
		#expand("balanced_cohorts/{seed}/{f}.tsv", seed=seeds, f=("metadata", "calls", "cohorts"))
	conda: "envs/r.yaml"
	params:
		size=config["balanced_cohort_size"],
		seeds=seeds,
		mem="4G",
		time="02:00:00"
	shell:
		"mapfile -t < {input.maxthresh} maxthresh && "
		#"for seed in {params.seeds}; do "
		"balanced_cohort_genomes.R "
		"--cohorts {input.cohorts} "
		"--threshold ${{maxthresh[0]}} "
		"--metadata {input.metadata} "
		"--size {params.size} "
		"--calls {input.calls} "
		"--seed {wildcards.seed} "
		"--outdir balanced_cohorts/{wildcards.seed}; "
		#"done"

use rule transpose_binary_markers as transpose_balanced with:
	input:
		"balanced_cohorts/{seed}/calls.tsv"
	output:
		"balanced_cohorts/{seed}/transposed_calls.tsv"

use rule run_feht as run_feht_balanced with:
	input:
		maxthresh=rules.find_maximum_negative_threshold.output,
		markers="balanced_cohorts/{seed}/transposed_calls.tsv",
		metadata="balanced_cohorts/{seed}/metadata.tsv",
		cohorts="balanced_cohorts/{seed}/cohorts.tsv"
	output:
		"balanced_cohorts/{seed}/feht_results.tsv"
	params:
		mem="2G",
		time="15:00"

use rule calculate_ntb as calculate_ntb_balanced with:
	input:
		"balanced_cohorts/{seed}/feht_results.tsv"
	output:
		"balanced_cohorts/{seed}/ntb.tsv"

rule symlink_balanced_ntb:
	input:
		expand("balanced_cohorts/{seed}/ntb.tsv", seed=seeds)
	output:
		expand("balanced_ntb/{seed}-ntb.tsv", seed=seeds)
	params:
		mem="2G",
		time="30:00"
	shell:
		"for seed_dir in balanced_cohorts/*/; do "
		"seed=$(basename $seed_dir); "
		"ln -sr $seed_dir/ntb.tsv balanced_ntb/${{seed}}-ntb.tsv; "
		"done"

rule summarize_balanaced:
	input:
		rules.symlink_balanced_ntb.output
	output:
		ranks="summary/ranks.tsv",
		ranks_summary="summary/ranks_summary.tsv",
		ntb="summary/marker-ntb.tsv",
		ntb_summary="summary/marker-ntb_summary.tsv"
	params:
		mem="4G",
		time="30:00"
	shell:
		"rank_markers.py --results-dir balanced_ntb/ --outdir summary/"

def get_pangenome_genes():
	genes=glob_wildcards("pangenome/feature_sequences/{gene}.nucleotide.fasta").gene
	return genes

rule generate_consensus_sequence:
	input:
		dummy=rules.pirate_pangenome.output,
		actual="pangenome/feature_sequences/{gene}.nucleotide.fasta"
	output:
		"consensus_sequences/{gene}.fasta"
	conda:
		"envs/emboss.yaml"
	threads:
		1
	params:
		mem="1G",
		time="45:00"
	shell:
		"cons -sequence {input.actual} -outseq {output} -name {wildcards.gene}"

rule generate_primer_configs:
	input:
		"consensus_sequences/{gene}.fasta"
	output:
		"primer3_configs/{size}/{gene}.txt"
	threads: 1
	params:
		mem="500M", time="10:00"
	shell:
		"format-primer3.py "
		"--input {input} "
		"--amplicon-size {wildcards.size} "
		"--primer-size 20 25 30 "
		"--output {output}"

rule generate_primer:
	input: "primer3_configs/{size}/{gene}.txt"
	output: "primers/{size}/{gene}.txt"
	conda: "envs/primer3.yaml"
	threads: 1
	params: mem="500M", time="10:00"
	shell: "primer3_core < {input} > {output}"

def aggregate_primers(wildcards):
	import pandas as pd
	nucleotide_fastas = list(Path("pangenome/feature_sequences/").glob("*.fasta"))
	markers = pd.read_csv("summary/ranks_summary.tsv", sep = "\t", index_col=0).index
	top100 = [m.split("_")[0] for m in markers][:100]
	genes = [m for m in top100 if Path(f"pangenome/feature_sequences/{m}.nucleotide.fasta").exists()]

	return expand("primers/{size}/{gene}.txt", size=primer_sizes, gene=genes)

rule format_primer_table:
	input: aggregate_primers #expand("primers/{size}/{gene}.txt", size=primer_sizes, gene=get_pangenome_genes())
	output:
		"results/primer_table.tsv"
	threads: 1
	params: mem="4G", time="30:00"
	shell: "tabulate-primers.py --input primers/ --output {output}"

rule select_top_primers:
	input: primers="results/primer_table.tsv", ranks="summary/ranks_summary.tsv"
	output: "results/top_primers.tsv"
	conda: "envs/r.yaml"
	threads: 1
	params: mem="1G", time="10:00", top_n=100
	shell: "select-ranked-primers.R {input.primers} {input.ranks} {params.top_n} {output}"

rule extract_top_pairs:
	input: "results/top_primers.tsv"
	output: "results/top_pairs.tsv"
	threads: 1
	params: mem="1G", time="10:00"
	shell: "xsv fmt -d '\t' {input} | xsv select 'LEFT_SEQUENCE,RIGHT_SEQUENCE' | xsv fmt -t '\t' | grep -v 'n' > {output}"

rule fastulate_top_pairs:
	input: "results/top_pairs.tsv"
	output: "results/top_pairs.fasta"
	threads: 1
	params: mem="1G", time="10:00"
	shell: "fastulate_primer_pairs.py --input {input} --output {output}"

rule blast_top_pairs_against_genomes:
	input: queries="results/top_pairs.fasta", genome="high_quality_genomes/{sample}.fasta"
	output: "blast_results/{sample}.tsv"
	conda: "envs/blast.yaml"
	threads: 1
	params: mem="1G", time="10:00"
	shell: "blastn -task blastn-short -subject {input.genome} -query {input.queries} -outfmt 6 > {output}"

def aggregate_genome_names():

	genome_fastas = Path("high_quality_genomes/").glob("*.fasta")
	genomes = [genome.stem for genome in genome_fastas]
	return genomes

rule concatenate_blast_results:
	input: expand("blast_results/{sample}.tsv", sample=aggregate_genome_names())
	output: "results/primer_blast.tsv"
	threads: 1
	params: mem="4G", time="30:00"
	shell: "cat {input} > {output}"

rule get_positive_primers:
	input: "results/primer_blast.tsv"
	output: "results/positive_hits.tsv"
	threads: 1
	params: mem="4G", time="10:00"
	shell: "xsv fmt -d '\t' {input} | xsv search -s length ^25$ | xsv search -s mismatches '[0-1]' > {output}"

rule generate_primer_pair_calls:
	input: "results/positive_hits.tsv"
	output: "results/primer_pair_calls.tsv"
	threads: 1
	params: mem="4G", time="10:00"
	shell: "generate-primer-pair-calls.R {input} {output}"

