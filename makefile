.PHONY: run build

run:
	@./build/src/RadixSort

build:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Release ../ -G Ninja; ninja;

debug:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Debug ../ -G Ninja; ninja;
