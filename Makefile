all:
	make --directory=src
	make --directory=exp

clean:
	make clean --directory=src
	make clean --directory=exp
