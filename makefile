RCODES := R/*.R

define test_file
    Rscript -e "devtools::test_active_file(file = './tests/testthat/test-$1.R')"
endef

t="sequence"

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

# useage
# make test_file t=file
test_file:
	$(call test_file,$t)
