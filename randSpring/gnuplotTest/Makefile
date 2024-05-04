graph_spring_hamiltonian: spring_hamiltonian
	./spring_hamiltonian

spring_hamiltonian:
	g++ spring_hamiltonian.cpp -o spring_hamiltonian -lboost_iostreams -lboost_system -lboost_filesystem

graph_spring_simple: spring_simple.dat
	gnuplot -p spring_simple.p

spring_simple.dat: spring_simple
	./spring_simple

spring_simple:
	g++ spring_simple.cpp -o spring_simple

clean:
	rm spring_simple spring_simple.dat -f
	rm spring_hamiltonian -f
