all :
	cd flat && $(MAKE)
	cd bw && $(MAKE)
	cd lval && $(MAKE)

clean :
	cd flat && $(MAKE) clean
	cd bw && $(MAKE) clean
	cd lval && $(MAKE) clean

.PHONY: all clean

