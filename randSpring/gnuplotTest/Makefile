graph_1D_spring_system: 1D_spring.mp4
	open 1D_spring.mp4

1D_spring.mp4: 1D_spring_system
	./1D_spring_system
	ffmpeg -pattern_type glob -i "./tmp/*.png" 1D_spring.mp4

1D_spring_system:
	g++ 1D_spring_system.cpp -o 1D_spring_system -lfftw3 -lm -lboost_iostreams -lboost_system -lboost_filesystem -llapacke -llapack -lblas

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
	rm 1D_spring_system -f
	rm *.mp4 ./tmp -rf
