main: sum test
	gcc *.o -o result
sum:
	gcc -c sum.c
test:
	gcc -c test.c
clean:
	rm -f *.o
	rm -f result
