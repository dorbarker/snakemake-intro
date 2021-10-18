from pathlib import Path
samples = [p.stem for p in Path("genomes").glob("*.fasta")]

rule all:
	input: "pangenome/PIRATE.gene_families.tsv"

rule annotate:
	input: "genomes/{sample}.fasta"
	output: "annotations/{sample}/{sample}.gff"
	threads: 8
	params: outdir="annotations/{sample}"
	shell:
		"prokka --force --cpus {threads} "
		"--prefix {wildcards.sample} --outdir {params.outdir} "
		"{input}"

rule symlink_gffs:
	input: "annotations/{sample}/{sample}.gff"
	output: "gffs/{sample}.gff"
	threads: 1
	shell: "ln -sr {input} {output}"

rule pangenome:
	input: expand("gffs/{sample}.gff", sample=samples)
	output: "pangenome/PIRATE.gene_families.tsv"
	threads: 8
	shell: "PIRATE --input gffs/ --output pangenome/ --nucl --threads {threads}"
