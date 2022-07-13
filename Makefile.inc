EXE=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
CPPC=g++ 
CPPC_FLAGS=-std=c++17 -Ofast 

all : 
	cp ../*.cpp .
	make exe
	time ./black_scholes_pde.exe 300 >> black_scholes_pde.out
	time ./burgers.exe 100 10000 >> burgers.out
	time ./libor.exe 10000 >> libor.out
	time ./black_scholes_mc.exe 10000000 >> black_scholes_mc.out
	./example.exe 3 1 
	dot -Tpdf dag.dot -o dag.pdf

exe : $(EXE)

%.exe : %.cpp
	$(CPPC) $(CPPC_FLAGS) $< -o $@ 

clean :
	rm -fr *.cpp *.exe *.dot *.pdf *.out

.PHONY: all exe clean 
