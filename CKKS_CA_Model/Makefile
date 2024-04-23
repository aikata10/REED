.PHONY = all clean

CC = g++

LINKERFLAG = -lm

SRCS := main.cpp
BINS := main

all: main

main: main.cpp
	@echo "Checking.."
	g++ main.cpp  -fopenmp -o main

clean:
	@echo "Cleaning up..."
	rm -rvf  main

