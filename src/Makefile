CC = gcc
C++FLAGS = -std=c++11 
GCOVFLAGS = -fprofile-arcs -ftest-coverage

OS := $(shell uname)

ifeq ($(OS), Darwin)
    FLAGSYSTEM = -lgtest
else
    FLAGSYSTEM = -pthread -lgtest -lstdc++ -lrt -lm
endif

all: test

s21_matrix_oop.a:
	$(CC) $(C++FLAGS) -c s21_matrix.cc s21_test_matrix.cc
	ar rcs s21_matrix_oop.a *.o
	ranlib s21_matrix_oop.a

test:  
	$(CC) $(C++FLAGS) s21_matrix.cc s21_test_matrix.cc -o s21_test $(FLAGSYSTEM)
	./s21_test


# gcov_report:    
# 	$(CC) $(GCOVFLAGS) $(CFLAGS) s21_matrix.cc s21_test_matrix.cc -o s21_test $(FLAGS)
# 	./s21_test  
# 	lcov -t "s21_test" -o test.info -c -d .
# 	genhtml -o report test.info
# 	rm -rf gcovreport gcovreport.info *.gcda *.gcno *.dSYM *info
# 	open report/index.html
rebuild: clean all

clean: 
	rm -rf *.o s21_test s21_report *.a a.out *.gcno *.gcda s21_test.dSYM *.info report 3DViewer dist s21_matrix
cf:
	clang-format --style=Google -n $(shell find . -name '*.c') $(shell find . -name '*.h') $(shell find . -name '*.cpp')
	clang-format --style=Google -i $(shell find . -name '*.c') $(shell find . -name '*.h') $(shell find . -name '*.cpp')
# install:
# 	mkdir 3DViewer
# 	cd 3DViewer && qmake ../My3DViewer && make
# open:
# 	./3DViewer/My3DViewer.app/Contents/MacOS/My3DViewer
# uninstall:
# 	rm -rf 3DViewer
# dvi:
# 	open ../README.md
# dist:
# 	mkdir dist
# 	cp -R My3DViewer Test_core Makefile dist
# 	cd dist && tar cvzf 3DViewer.tgz *
# 	cd dist && rm -rf *.c *.cpp *.h *.ui calc* My3DViewer Makefile Test_core

