
rule all:
	input: expand("presentations/snakemake-intro-lesson-{lesson}.{ext}", lesson=[1,2,3], ext=["pdf", "html"])

rule render_html:
	input:
		"lessons/{filename}.Rmd"
	output:
		"presentations/{filename}.html"
	shell:
		"R --quiet -e 'rmarkdown::render(\"{input}\", output_dir=\"presentations\")'"

rule render_pdf:
	input:
		"lessons/{filename}.Rmd"
	output:
		"presentations/{filename}.pdf"
	priority: 1
	shell:
		"R --quiet -e 'pagedown::chrome_print(\"{input}\", output=\"{output}\")'"
