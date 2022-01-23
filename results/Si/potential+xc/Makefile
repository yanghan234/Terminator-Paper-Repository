src = real.cpp
target = real.x
headers = -I/usr/local/include
libflag = -L/usr/local/lib -lgsl -lgslcblas

${target} : ${src}
	g++  -std=c++11 -o $@ $^ ${headers} ${libflag}
