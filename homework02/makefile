CFLAGS = -fopenmp -O3
CC = g++-9

all: val_test01_solved val_test02_solved omp_bug2 omp_bug3 omp_bug4 omp_bug5 omp_bug6 gs ja

val_test01_solved: val_test01_solved.cpp
	${CC} -o val_test01_solved ${CFLAGS} val_test01_solved.cpp

val_test02_solved: val_test02_solved.cpp
	${CC} -o val_test02_solved ${CFLAGS} val_test02_solved.cpp

omp_bug2: omp_bug2_solved.c 
	${CC} -o omp_solved2 ${CFLAGS} omp_bug2_solved.c 

omp_bug3: omp_bug3_solved.c 
	${CC} -o omp_solved3 ${CFLAGS} omp_bug3_solved.c 

omp_bug4: omp_bug4_solved.c 
	${CC} -o omp_solved4 ${CFLAGS} omp_bug4_solved.c 

omp_bug5: omp_bug5_solved.c 
	${CC} -o omp_solved5 ${CFLAGS} omp_bug5_solved.c 

omp_bug6: omp_bug6_solved.c 
	${CC} -o omp_solved6 ${CFLAGS} omp_bug6_solved.c 

gs: gs2D-omp.cpp residual.h
	${CC} -o gs2D-omp ${CFLAGS} residual.h gs2D-omp.cpp

ja: jacobi2D-omp.cpp residual.h
	${CC} -o jacobi2D-omp ${CFLAGS} residual.h jacobi2D-omp.cpp