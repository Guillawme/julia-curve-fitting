PROJECT = julia-curve-fitting

FILES = \
	Project.toml \
	Manifest.toml \
	build-sysimg.jl \
	makefile

.PHONY: all 

all: sysimg.dylib

sysimg.dylib: $(FILES)
	julia build-sysimg.jl
