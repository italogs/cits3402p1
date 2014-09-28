run: compilado
	./outputModified intel
compilado:
	gcc -fopenmp -std=c99 myCode.c -o outputModified
