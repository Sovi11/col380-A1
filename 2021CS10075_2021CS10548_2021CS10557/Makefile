all:
	@echo "Pthreads :"
	@make pthreads_array
	@echo ""
	@echo "OpenMP :" 
	@make omp_array
	@echo ""
	@echo "Sequential :"
	@make sequential_array
	@echo ""

all_optimized:
	@echo "Pthreads :"
	@make pthreads_array_optimized
	@echo ""
	@echo "OpenMP :" 
	@make omp_array_optimized
	@echo ""
	@echo "Sequential :"
	@make sequential_array_optimized
	@echo ""

pthreads_array : 
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 pth_arrays.cpp -o pth_arrays.o -lpthread
	@./pth_arrays.o $(DIM) $(THREADS) $(CHCK)

omp_array : 
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 omp_array.cpp -o omp_array.o -fopenmp
	@./omp_array.o $(DIM) $(THREADS) $(CHCK)

sequential_array :
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 sequential_array.cpp -o sequential_array.o
	@./sequential_array.o $(DIM) $(CHCK)

pthreads_array_optimized : 
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 -O3 pth_arrays.cpp -o pth_arrays.o -lpthread
	@./pth_arrays.o $(DIM) $(THREADS) $(CHCK)

omp_array_optimized : 
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 -O3  omp_array.cpp -o omp_array.o -fopenmp
	@./omp_array.o $(DIM) $(THREADS) $(CHCK)

sequential_array_optimized :
	@python3 generate.py $(DIM) > input.txt 
	@g++-13 -O3 sequential_array.cpp -o sequential_array.o
	@./sequential_array.o $(DIM) $(CHCK)

clean : 
	@rm -f *.o input.txt
