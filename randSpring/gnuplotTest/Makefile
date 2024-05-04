graph_symplectic_spring: symplectic_spring
	./symplectic_spring

symplectic_spring:
	g++ symplectic_spring.cpp -o symplectic_spring -lboost_iostreams -lboost_system -lboost_filesystem

graph_euler_spring: euler_spring.dat
	gnuplot -p euler_spring.p

euler_spring.dat: euler_spring
	./euler_spring

euler_spring:
	g++ euler_spring.cpp -o euler_spring

clean:
	rm euler_spring euler_spring.dat -f
	rm symplectic_spring -f
