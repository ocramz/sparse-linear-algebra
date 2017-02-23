.DEFAULT_GOAL := help

# number of CPU threads
NTHREADS=2

help:
	@echo "  Use \`make <target>\` where <target> is one of"
	@echo "  clean          remove .stack-work/"
	@echo "  build          build and link (`stack build`)"
	@echo "  build-profile  (`stack build --profile`)"
	@echo "  test           (`stack test`)"
	@echo "  test-profile   (`stack test --profile`)"


clean:
	rm -rf .stack-work/

build:
	stack build

test:
	stack test

build-threaded:
	stack build --ghc-options "-threaded -O2"


build-profile:
	stack build --profile
	# stack build --executable-profiling --library-profiling --ghc-options="-fprof-auto -rtsopts"

test-profile:
	# stack test --profile
        # stack exec $(stack path --dist-dir)/build/spec -- +RTS -h
	stack test --ghc-options "+RTS -h -RTS"

test-threaded:
	stack test --ghc-options "+RTS -N${NTHREADS} -s -RTS"

all:
	make build
	make test


rebuild:
	make clean
	make all

profile:
	make clean
	make build-profile
	make test-profile

threaded:
	make build-threaded
	make test-threaded
