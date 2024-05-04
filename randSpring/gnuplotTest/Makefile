graph_spring_simple: spring_simple.dat
	gnuplot -p spring_simple.p

spring_simple.dat: spring_simple
	./spring_simple

spring_simple:
	g++ spring_simple.cpp -o spring_simple

clean:
	rm spring_simple -f
