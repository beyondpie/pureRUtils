RCODES := R/*.R


install:
	Rscript -e "remotes::install_github('beyondpie/pureRUtils')"

all: README.md doc

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

doc: $(RCODES)
	Rscript -e "roxygen2::roxygenize()"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test(pkg = '.')"
