from pathlib import Path
samples = [p.stem for p in Path("genomes").glob("*.fasta")]

gbk_file = config["proteins"]

localrules: all, symlink_gffs

rule all:
	input: "core_loci.txt"

rule annotate:
	input: "genomes/{sample}.fasta"
	output: "annotations/{sample}/{sample}.gff"
	conda: "envs/prokka.yaml"
	threads: 8
	params: mem="12G", time="20:00"
	shell:
		"prokka --force --cpus {threads} "
		"--prefix {wildcards.sample} --outdir annotations/{wildcards.sample} "
		"--proteins {gbk_file} "
		"{input}"

rule symlink_gffs:
	input: "annotations/{sample}/{sample}.gff"
	output: "gffs/{sample}.gff"
	threads: 1
	shell: "ln -sr {input} {output}"

rule pangenome:
	input: expand("gffs/{sample}.gff", sample=samples)
	output: "pangenome/PIRATE.gene_families.tsv"
	conda: "envs/pirate.yaml"
	threads: 8
	params: mem="12G", time="30:00"  # In real life, pirate takes MUCH MORE than this
	shell: "PIRATE --input gffs/ --output pangenome/ --nucl --threads {threads}"

rule get_core_loci:
	input: "pangenome/PIRATE.gene_families.tsv"
	output: "core_loci.txt"
	params: mem="4G", time="10:00"
	run:
		import pandas as pd
		pangenome = pd.read_csv(input[0], sep="\t")
		core_genes = pangenome["number_genomes"] == 12
		selected_loci_names = pangenome["gene_family"].loc[core_genes]
		selected_loci_names.to_csv(output[0], header=False)

