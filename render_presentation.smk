rule all:
	input: "snakemake-intro.html", "snakemake-intro.pdf"

rule render_html:
	input:
		"{filename}.Rmd"
	output:
		"{filename}.html"
	shell:
		"R --quiet -e 'rmarkdown::render(\"{input}\")'"

rule render_pdf:
	input:
		"{filename}.Rmd"
	output:
		"{filename}.pdf"
	shell:
		"R --quiet -e 'pagedown::chrome_print(\"{input}\")'"
