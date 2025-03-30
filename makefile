.PHONY: run build debug

run:
	@./build/src/RadixSort 5 3

build:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Release ../ -G Ninja; ninja;

debug:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Debug ../ -G Ninja; ninja;
