AUTHOR = ocramz
REPO = petsc-hs

.DEFAULT_GOAL := help



help:
	@echo "Use \`make <target>\` where <target> is one of"
	@echo "  clean          remove .stack-work/"
	@echo "  build          build and link (`stack build`)"
	@echo "  build-profile  (`stack build --profile`)"
	@echo "  test           (`stack test`)"
	@echo "  test-profile   (`stack test --profile`)"


clean:
	rm -rf .stack-work/

build:
	stack build

build-profile:
	stack build --profile

test-profile:
	stack test --profile


all:
	make clean
	make build
	make test

profile:
	make clean
	make build-profile
	make test-profile
