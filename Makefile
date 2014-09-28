run: compilado
	./outputMod intel
compilado:
	gcc -fopenmp -std=c99 myCode.c -o outputMod