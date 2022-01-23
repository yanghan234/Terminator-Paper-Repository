CXX = g++ -std=c++11
INCS = -I/usr/include
LIBS = -L/usr/local/lib -lgsl -lgslcblas
SRC = ldos.cpp

ldos.x : ${SRC}
	${CXX} -o $@ ${LIBS} ${INCS}  $^

clean :
	rm *.x
