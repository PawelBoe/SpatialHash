SOURCES = $(wildcard ./*.cpp)

all: benchmark test

benchmark:
	g++ -O3 -std=c++17 -Wall -I .. benchmark.cpp -o hash_benchmark

test:
	g++ -O3 -std=c++17 -Wall -I .. test.cpp -o spatial_hash_test

clean:
	rm spatial_hash_test hash_benchmark
