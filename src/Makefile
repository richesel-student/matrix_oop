CC = gcc -g
CFLAGS = -Wall -Werror -Wextra -std=c++17 -pedantic
CPPFLAGS = -lstdc++
GCOV = --coverage

OS = $(shell uname)
ifeq ($(OS), Linux)
	TEST = -lgtest -lm -lsubunit -lpthread
else
	TEST = -lgtest -lm -lpthread
endif

all: s21_matrix_oop.a

clean:
	rm -rf *.o *.a test *.gcno *.gcda *.info report

test: s21_matrix_oop.a test_s21_matrix_oop.cpp 
	$(CC) $(CFLAGS) -c test_s21_matrix_oop.cpp 
	$(CC) $(CFLAGS) $(GCOV) *.o $(CPPFLAGS) $< $(TEST) -o $@
	./$@

s21_matrix_oop.a: s21_matrix_oop.cpp  s21_matrix_oop.h
	$(CC) $(CFLAGS) $(GCOV) -c s21_matrix_oop.cpp 
	ar rc $@ *.o
	ranlib $@

gcov_report: test
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

check_style:
	clang-format -n -style=google *.h *.cpp cases/*.cpp

correct_style:
	clang-format -i -style=google *.h *.cpp cases/*.cpp

cppcheck:
	cppcheck --enable=all --inconclusive --std=c++17 --language=c++ --suppress=missingIncludeSystem *.h *.cpp cases/*.cpp

valgrind:
	make test
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -s --verbose -q --show-reachable=yes --num-callers=20 --trace-children=yes --track-fds=yes ./test
 #valgrind --leak-check=full ./test
	
leaks:
	make test
	leaks -quiet -atExit -- ./test

.PHONY: all clean test s21_matrix_oop.a gcov_report check_style correct_style cppcheck valgrind leaks