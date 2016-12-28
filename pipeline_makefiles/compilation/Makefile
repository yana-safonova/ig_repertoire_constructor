include ../Makefile_vars

all: compile

compile:
	cd $(IGREC) && ./prepare_cfg
	$(MAKE) -C $(IGREC) -j7

clean:
	$(MAKE) -C $(IGREC) clean
