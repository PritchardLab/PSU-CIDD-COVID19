.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

all:
	$(MAKE) -C cpp-*
	mkdir -p outputs
	cd cpp-* && ./basic_traj_graph.r
	mv cpp-*/basic_traj.png outputs/

clean:
	$(MAKE) -C cpp-v4-6e-severe-classes clean


