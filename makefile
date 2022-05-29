# NOTE
# On my MacOS, I installed R with brew, not with conda.
# But conda base will change the PATH, which cause errors.
# So need to update the PATH by removing conda-related variables.

RCODES := R/*.R

define test_file
    Rscript -e "devtools::test_active_file(file = './tests/testthat/test-$1.R')"
endef

t="sequence"

install:
	Rscript -e "remotes::install_github('beyondpie/pureRUtils', upgrade = 'never')"

all: README.md doc

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

doc: $(RCODES)
	Rscript -e "roxygen2::roxygenize()"

check:
	Rscript -e "devtools::check(error_on = 'never')"

test:
	Rscript -e "devtools::test(pkg = '.')"

# useage
# make test_file t=file
test_file:
	$(call test_file,$t)
