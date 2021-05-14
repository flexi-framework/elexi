# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_userguide.md 010_gettingstarted.md 020_codeoptions.md 021_parameter.md  030_workflow.md 040_toolsoverview.md
SRC +=  050_tutorials.md tutorials/*/*.md 070_unittest.md 080_installationguidelines.md
ifneq ($(BUILDALL),)
  SRC += missing.md
endif
SRC += references.md

# Pandoc changed citeproc in 2.8
PANDOCVERSION := $(shell pandoc --version | grep ^pandoc | sed 's/^.* //g')
PANDOCMAJOR   := $(shell echo $(PANDOCVERSION) | cut -f1 -d.)
PANDOCMINOR   := $(shell echo $(PANDOCVERSION) | cut -f2 -d.)
ifeq ($(shell expr $(PANDOCMAJOR) \= 2), 0)
  $(error No or unkown version of pandoc. Please use pandoc v2.x to build this documentation.)
endif

ifeq ($(shell expr $(PANDOCMINOR) \>= 8), 1)
  PANDOCFILTER=--citeproc
else
  PANDOCFILTER=--filter pandoc-citeproc
endif

# Update template to current version if possible
TEMPLATE_INFO := $(shell pandoc -D latex > /dev/null; echo $$?)
ifeq ($(shell expr $(TEMPLATE_INFO) \= 0), 1)
  $(info Found pandoc latex template, using system version)
  $(shell pandoc -D latex > template.tex)
else
  $(info pandoc latex template not found on system, using git version)
endif

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
	pandoc $(SRC)   -o userguide.pdf -s -N $(PANDOCFILTER) --from markdown+implicit_figures --bibliography userguide.bib --listings --template=template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N $(PANDOCFILTER) --from markdown+implicit_figures --bibliography userguide.bib --listings --template=template.tex

.PHONY: userguide

