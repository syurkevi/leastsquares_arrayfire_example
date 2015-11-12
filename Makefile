INC =  -I /usr/local/include
LIBS= -lafcpu
LIB_PATHS= -L /usr/local/lib/
FLAGS =  -m64 -Wall -O3 -DNDEBUG -std=c++11

all :
	g++ -o lsq main.cpp ${FLAGS} ${INC} ${LIB_PATHS} ${LIBS}  
