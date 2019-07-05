########### Gnu:
#COMPILER = g++
#FLAGS= -static -Ofast -fipa-pta -fopenmp
#LIB = -lm -lgsl -lgslcblas -lgomp -pthread

############ Intel:

#COMPILER = icpc
#FLAGS=     -ipo -O3  -no-prec-div -xHost -simd -qopenmp -fp-model fast=2  -g
#LIB = -mkl -lgsl 

############ Intel 2 :

#COMPILER = icpc
#FLAGS= -ipo -O3  -no-prec-div -xHost -axAVX,CORE-AVX-I -simd -qopenmp -fp-model fast=2 -static
#LIB = -mkl -lgsl 

############ Intel : no AVX2

COMPILER = icpc
FLAGS= -ipo -O3  -xHost -simd -fp-model fast=2 -static
LIB = -mkl -lgsl 

############# Pgi:
#FLAGS =  -fast -Mipa=fast,inline -O3
CPPS = $(wildcard src/*.cpp)
HEADER = $(wildcard include/*.h)
OBJS  =$(patsubst src/%.cpp,.build/%.o,${CPPS})

castle: ${OBJS}
	@${COMPILER}  ${FLAGS}  ${OBJS} ${LIB} -o castle
${OBJS}: .build/%.o: src/%.cpp | .build
	${COMPILER}  ${FLAGS} -c $< -o $@
${OBJS}: ${HEADER}
.build:
	mkdir .build

clean:
	@rm -f uno.out
	@rm -rf .build
	@rm -f remove*
