.PHONY: run build

run:
	@./build/src/RadixSort

build:
	@mkdir -p build/; cd build/; cmake ../ -G Ninja; ninja;