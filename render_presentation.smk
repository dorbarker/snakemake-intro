
rule all:
	input: expand("lessons/snakemake-intro-lesson-{lesson}.{ext}", lesson=[1,2,3], ext=["pdf", "html"])

rule render_html:
	input:
		"lessons/{filename}.Rmd"
	output:
		"lessons/{filename}.html"
	shell:
		"R --quiet -e 'rmarkdown::render(\"{input}\")'"

rule render_pdf:
	input:
		"lessons/{filename}.Rmd"
	output:
		"lessons/{filename}.pdf"
	priority: 1
	shell:
		"R --quiet -e 'pagedown::chrome_print(\"{input}\", output=\"{output\"})'"
