LAFLAG=-llapacke -llapack -lblas
GFLAG=-lboost_iostreams -lboost_system -lboost_filesystem

play: nonlinear.mp4
	open nonlinear.mp4

run: nonlinear_spring_system
	./nonlinear_spring_system

Test_linalg: testLinalg
	./testLinalg

Test_graphSystem: testGraphSystem
	./testGraphSystem

nonlinear.mp4: nonlinear_spring_system
	./nonlinear_spring_system
	ffmpeg -pattern_type glob -i "./tmp/*.png" nonlinear.mp4

nonlinear_spring_system: linalg.o graph_system.o
	g++ nonlinear_spring_system.cpp linalg.o graph_system.o -o nonlinear_spring_system $(LAFLAG) $(GFLAG)

testLinalg: linalg.o
	g++ testLinalg.cpp linalg.o -o testLinalg $(LAFLAG)

testGraphSystem: graph_system.o linalg.o
	g++ testGraphSystem.cpp linalg.o graph_system.o -o testGraphSystem $(LAFLAG) $(GFLAG)

graph_system.o:
	g++ -c graph_system.cpp -o graph_system.o

linalg.o:
	g++ -c linalg.cpp -o linalg.o

clean:
	rm linalg.o testLinalg -f
	rm graph_system.o testGraphSystem -f
	rm linalg.o testLinalg -f
	rm nonlinear_spring_system -f
	rm *.mp4 ./tmp -rf
