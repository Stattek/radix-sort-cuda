.PHONY: run build

run:
	@./build/RadixSort

build:
	@mkdir -p build/; cd build/; cmake ../ -G Ninja; ninja;