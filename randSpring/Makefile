run: testLinalg
	./testLinalg

testLinalg: linalg.o
	g++ testLinalg.cpp linalg.o -o testLinalg -llapacke -llapack -lblas

linalg.o:
	g++ -c linalg.cpp -o linalg.o

clean:
	rm linalg.o testLinalg -f
