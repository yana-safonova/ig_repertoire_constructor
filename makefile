rig:
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dig:
	$(MAKE) -C build/debug/ig_repertoire_constructor ig_repertoire_constructor

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all

