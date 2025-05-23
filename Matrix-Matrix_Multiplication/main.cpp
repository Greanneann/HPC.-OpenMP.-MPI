#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "include/Stopwatch.hpp"
#include <omp.h>


#define CACHE_CONST_BLOCKING_SIZE 16


/**
 * Different implementations
 *
 * Add your own test case here!
 */
enum
{
	MATRIX_SUM_ROWWISE = 10,
	MATRIX_SUM_COLWISE = 11,

	MATRIX_MATRIX_MUL_SIMPLE_IJK = 20,
	MATRIX_MATRIX_MUL_SIMPLE_JIK = 21,
	MATRIX_MATRIX_MUL_SIMPLE_IKJ = 22,
	MATRIX_MATRIX_MUL_SIMPLE_JKI = 23,
	MATRIX_MATRIX_MUL_SIMPLE_KIJ = 24,
	MATRIX_MATRIX_MUL_SIMPLE_KJI = 25,

	MATRIX_MATRIX_MUL_RESTRICTED_IKJ = 31,
	MATRIX_MATRIX_MUL_VAR_BLOCKED_IKJ = 32,
	MATRIX_MATRIX_MUL_BLOCKED_IKJ = 33,
	MATRIX_MATRIX_MUL_SIMD_IKJ = 34,
	MATRIX_MATRIX_MUL_OPENMP_IKJ = 35
};


/**
 * Output how to use the program
 */
void print_program_usage(int argc, char *argv[])
{
	std::cout << std::endl;
	std::cout << "Program usage:" << std::endl;
	std::cout << std::endl;
	std::cout << "  " << argv[0] << " [mat-mat-mul variant (int)] [N probem size (int) >= 1] [cache block size (int) >= 1] " << std::endl;
	std::cout << std::endl;
}


/**
 * Output the matrix content
 */
void print_matrix(std::size_t N, const double *C)
{
	for (std::size_t i = 0; i < N; i++)
	{
		for (std::size_t j = 0; j < N; j++)
		{
			std::cout << C[i*N+j];
			if (j < N-1)
				std::cout << "\t";
		}
		std::cout << std::endl;
	}
}


double kernel__matrix_sum_rowwise(std::size_t N, const double *i_A)
{
	/*
	 * STUDENT ASSIGNMENT
	 */
	double acc = 0;
	// A[i][j] is equivalent to A[i*N + j], where i is a row, and j is a column 
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t j = 0; j < N; j++){
			acc += i_A[i * N + j];
		}
	}
	return acc;
}


double kernel__matrix_sum_colwise(std::size_t N, const double *i_A)
{
	/*
	 * STUDENT ASSIGNMENT
	 */
	// A[i][j] is equivalent to A[i*N + j], where i is a row, and j is a column 
	double acc = 0;
	for(std::size_t j = 0; j < N; j++){
		for(std::size_t i = 0; i < N; i++){
			acc += i_A[i * N + j]; 
		}
	}
	return acc;
}


/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_ijk(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t j = 0; j < N; j++){
			double sum = 0; 
			for(std::size_t k = 0; k < N; k++){
				sum += i_A[i*N + k]*i_B[k*N + j];
			}
			o_C[i*N + j] += sum; 
		}
	}
}


/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_jik(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t j = 0; j < N; j++){
		for(std::size_t i = 0; i < N; i++){
			double sum = 0; 
			for(std::size_t k = 0; k < N; k++){
				sum += i_A[i*N + k]*i_B[k*N + j];
			}
			o_C[i*N + j] += sum; 
		}
	}
}



/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_jki(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t j = 0; j < N; j++){
		for(std::size_t k = 0; k < N; k++){
			for(std::size_t i = 0; i < N; i++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}


/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_kij(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t k = 0; k < N; k++){
		for(std::size_t i = 0; i < N; i++){
			for(std::size_t j = 0; j < N; j++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}


/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_kji(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t k = 0; k < N; k++){
		for(std::size_t j = 0; j < N; j++){
			for(std::size_t i = 0; i < N; i++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}

/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simple_ikj(std::size_t N, const double *i_A, const double *i_B, double *o_C)
{
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t k = 0; k < N; k++){
			// #pragma omp simd aligned(i_A, i_B, o_C : 512)
			for(std::size_t j = 0; j < N; j++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}

/**
 * Run restricted matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_restricted_ikj(std::size_t N, const double * __restrict__ i_A, const double * __restrict__ i_B, double * __restrict__ o_C)
{
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t k = 0; k < N; k++){
			// #pragma omp simd aligned(i_A, i_B, o_C : 512)
			for(std::size_t j = 0; j < N; j++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}



/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_var_blocked_ikj(std::size_t N, const double * __restrict__ i_A, const double * __restrict__ i_B, double * __restrict__ o_C, std::size_t cache_blocking_size)
{
    for (std::size_t ii = 0; ii < N; ii += cache_blocking_size) { // start at particular place in original matrix for i index
        for (std::size_t kk = 0; kk < N; kk += cache_blocking_size) { // start at particular place in original matrix for k index
            for (std::size_t jj = 0; jj < N; jj += cache_blocking_size) { // start at particular place in original matrix for j index

                // loops inside the particular block (ii,kk) and (kk, jj)
                for (std::size_t i = ii; i < std::min(ii + cache_blocking_size, N); i++) { // i < std::min(ii + cache_blocking_size, N) not to go outside the matri
                    for (std::size_t k = kk; k < std::min(kk + cache_blocking_size, N); k++) {
						// #pragma omp simd aligned(i_A, i_B, o_C : 512)
                        for (std::size_t j = jj; j < std::min(jj + cache_blocking_size, N); j++) {
                            o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
                        }
                    }
                }
            }
        }
    }
}



/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_blocked_ikj(std::size_t N, const double * __restrict__ i_A, const double * __restrict__ i_B, double * __restrict__ o_C)
{
    for (std::size_t ii = 0; ii < N; ii += CACHE_CONST_BLOCKING_SIZE) { // start at particular place in original matrix for i index
        for (std::size_t kk = 0; kk < N; kk += CACHE_CONST_BLOCKING_SIZE) { // start at particular place in original matrix for k index
            for (std::size_t jj = 0; jj < N; jj += CACHE_CONST_BLOCKING_SIZE) { // start at particular place in original matrix for j index

                // loops inside the particular block (ii,kk) and (kk, jj)
                for (std::size_t i = ii; i < std::min(ii + CACHE_CONST_BLOCKING_SIZE, N); i++) { // i < std::min(ii + cache_blocking_size, N) not to go outside the matri
                    for (std::size_t k = kk; k < std::min(kk + CACHE_CONST_BLOCKING_SIZE, N); k++) {
						// #pragma omp simd aligned(i_A, i_B, o_C : 512)
                        for (std::size_t j = jj; j < std::min(jj + CACHE_CONST_BLOCKING_SIZE, N); j++) {
                            o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
                        }
                    }
                }
            }
        }
    }
}


/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_simd_ikj(std::size_t N, const double * __restrict__ i_A, const double * __restrict__ i_B, double * __restrict__ o_C)
{
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t k = 0; k < N; k++){
			#pragma omp simd aligned(i_A, i_B, o_C : 512)
			for(std::size_t j = 0; j < N; j++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}

/**
 * Run simple matrix-matrix multiplication
 */
void kernel__matrix_matrix_mul_openmp_ikj(std::size_t N, const double * __restrict__ i_A, const double * __restrict__ i_B, double * __restrict__ o_C)
{
	#pragma omp parallel for collapse(2)
	for(std::size_t i = 0; i < N; i++){
		for(std::size_t k = 0; k < N; k++){
			for(std::size_t j = 0; j < N; j++){
				o_C[i*N + j] += i_A[i*N + k]*i_B[k*N + j];
			}
		}
	}
}

/**
 * Quick validation of the matrix C
 */
void validate_matrix_C(std::size_t N, double *C)
{
	for (std::size_t i = 0; i < N; i++)
	{
		for (std::size_t j = 0; j < N; j++)
		{
			double error = -1;
			double expected_value;
			if (i == j)
			{
				error = std::abs(C[i*N+j]-1);
				expected_value = 1.0;
			}
			else
			{
				error = std::abs(C[i*N+j]);
				expected_value = 0.0;
			}

			if (error > 1e-10*std::sqrt((double)N))
			{
				std::cerr << "******************************************************************" << std::endl;
				std::cerr << "C matrix has invalid entry" << std::endl;
				std::cerr << "******************************************************************" << std::endl;
				std::cerr << " + row: " << i << std::endl;
				std::cerr << " + col: " << j << std::endl;
				std::cerr << " + value: " << C[i*N+j] << std::endl;
				std::cerr << " + expected value: " << expected_value << std::endl;
				std::cerr << "******************************************************************" << std::endl;
				exit(1);
			}
		}
	}
}


/**
 * We initialize A by a trigonometric function
 */
void matrix_setup_A(std::size_t N, double *M)
{
	for (std::size_t i = 0; i < N; i++)
	{
		for (std::size_t j = 0; j < N; j++)
		{
            M[i * N + j] = cos((j + 0.5) * (i + 0.5) * M_PI / N); 
		}
	}
}


/**
 * We initialize B by a trigonometric function
 */
void matrix_setup_B(std::size_t N, double *M)
{
	/*
	 * STUDENT ASSIGNMENT
	 */
	for (std::size_t i = 0; i < N; i++)
	{
		for (std::size_t j = 0; j < N; j++)
		{
			M[i * N + j] = cos((j + 0.5) * (i + 0.5) * M_PI / N) * 2 / N;
		}
	}
}


/**
 * We initialize C by simply setting everything to -1
 * This data will be overwritten later on (hopefully).
 */
void matrix_zero_C(std::size_t N, double *M)
{
	for (std::size_t i = 0; i < N; i++)
	{
		for (std::size_t j = 0; j < N; j++)
		{
			M[j*N + i] = 0; // j is a number of row; i is a number of column
		}
	}
}


/**
 * Perform some dummy operations in memory to make sure
 * that the entire cache is not filled with any matrix
 * data
 */
void flush_cache()
{
	std::size_t N = 1024*1024*256/sizeof(double);	// 256 MB
	volatile double *buffer = new double[N];

	for (int i = 0; i < N; i++)
		buffer[i] = i;

	for (int i = 1; i < N; i++)
		buffer[i] = buffer[i-1];

	delete [] buffer;
}


double run_benchmark(int variant_id, long N, double *A, double *B, double*C, long cache_blocking_size)
{
	double retscalar = -1;
	switch (variant_id)
	{
		default:
			std::cerr << "Kernel not implemented" << std::endl;
			exit(-1);
			break;

		case MATRIX_SUM_COLWISE:
			retscalar = kernel__matrix_sum_colwise(N, A);
			break;

		case MATRIX_SUM_ROWWISE:
			retscalar = kernel__matrix_sum_rowwise(N, A);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IJK:
			kernel__matrix_matrix_mul_simple_ijk(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_JIK:
			kernel__matrix_matrix_mul_simple_jik(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IKJ:
			kernel__matrix_matrix_mul_simple_ikj(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_JKI:
			kernel__matrix_matrix_mul_simple_jki(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_KIJ:
			kernel__matrix_matrix_mul_simple_kij(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_KJI:
			kernel__matrix_matrix_mul_simple_kji(N, A, B, C);
			break;


		case MATRIX_MATRIX_MUL_RESTRICTED_IKJ:
			kernel__matrix_matrix_mul_restricted_ikj(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_VAR_BLOCKED_IKJ:
			kernel__matrix_matrix_mul_var_blocked_ikj(N, A, B, C, cache_blocking_size);
			break;


		case MATRIX_MATRIX_MUL_BLOCKED_IKJ:
			kernel__matrix_matrix_mul_blocked_ikj(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_SIMD_IKJ:
			kernel__matrix_matrix_mul_simd_ikj(N, A, B, C);
			break;

		case MATRIX_MATRIX_MUL_OPENMP_IKJ:
			kernel__matrix_matrix_mul_openmp_ikj(N, A, B, C);
			break;
	}

	return retscalar;
}


int main(int argc, char *argv[])
{
	std::cout << std::setprecision(10);
	std::cerr << std::setprecision(10);

	long N = -1;
	int variant_id = 1;
	long cache_blocking_size = CACHE_CONST_BLOCKING_SIZE;

	/**
	 * Variant ID
	 */
	if (argc >= 2)
	{
		variant_id = atoi(argv[1]);
	}

	/**
	 * Problem size
	 */
	if (argc >= 3)
	{
		N = atoi(argv[2]);
	}

	/**
	 * Bogus parameter which can be used for different things (e.g. blocking)
	 */
	if (argc >= 4)
	{
		cache_blocking_size = atoi(argv[3]);
	}

	if (N <= 0)
	{
		print_program_usage(argc, argv);
		return -1;
	}

	std::cout << " + variant_id: " << variant_id << std::endl;
	std::cout << " + N: " << N << std::endl;
	std::cout << " + cache_blocking_size: " << cache_blocking_size << std::endl;
	std::cout << " + size_per_matrix: " << (N*N*sizeof(double)) << std::endl;
	std::cout << " ++ k_size_per_matrix: " << (N*N*sizeof(double)*1e-3) << std::endl;
	std::cout << " ++ m_size_per_matrix: " << (N*N*sizeof(double)*1e-6) << std::endl;
	std::cout << " ++ g_size_per_matrix: " << (N*N*sizeof(double)*1e-9) << std::endl;
	std::cout << " ++ t_size_per_matrix: " << (N*N*sizeof(double)*1e-12) << std::endl;


    // Alignment for vectorized operations
    const std::size_t alignment = 512; 
    const std::size_t size = N * N * sizeof(double);

    double *A, *B, *C;
	
    // Allocate aligned memory for arrays A, B, and C
    if (posix_memalign((void**)&A, alignment, size) != 0 ||
        posix_memalign((void**)&B, alignment, size) != 0 ||
        posix_memalign((void**)&C, alignment, size) != 0) {
        std::cerr << "Error: Memory alignment failed" << std::endl;
        return 1;
    }

    // Alternative memory allocation without alignment (commented out)
    // Uncomment these lines if alignment is not required:

    // Allocate memory dynamically for arrays A, B, and C
    // A = new double[N*N];
    // B = new double[N*N];
    // C = new double[N*N];

	/**
	 * Setup particular benchmark
	 */
	const char *kernel_str = 0;

	switch (variant_id)
	{
		default:
			std::cerr << "***" << std::endl;
			std::cerr << "Setup not implemented" << std::endl;
			std::cerr << "***" << std::endl;
			exit(-1);
			break;

		case MATRIX_SUM_ROWWISE:
			kernel_str = "matrix_sum_rowwise";
			break;

		case MATRIX_SUM_COLWISE:
			kernel_str = "matrix_sum_colwise";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IJK:
			kernel_str = "mm_mul_simple_ijk";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_JIK:
			kernel_str = "mm_mul_simple_jik";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IKJ:
			kernel_str = "mm_mul_simple_ikj";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_JKI:
			kernel_str = "mm_mul_simple_jki";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_KIJ:
			kernel_str = "mm_mul_simple_kij";
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_KJI:
			kernel_str = "mm_mul_simple_kji";
			break;


		case MATRIX_MATRIX_MUL_RESTRICTED_IKJ:
			kernel_str = "mm_mul_restricted_ikj";
			break;


		case MATRIX_MATRIX_MUL_VAR_BLOCKED_IKJ:
			kernel_str = "mm_mul_var_blocked_ikj";

			if (N < cache_blocking_size)
			{
				std::cerr << "Cache blocking size must be at least N" << std::endl;
				exit(1);
			}
			break;

		case MATRIX_MATRIX_MUL_BLOCKED_IKJ:
			kernel_str = "mm_mul_blocked_ikj";

			if (N < CACHE_CONST_BLOCKING_SIZE)
			{
				std::cerr << "Precompiler cache blocking size must be at least N" << std::endl;
				exit(1);
			}
			break;

		case MATRIX_MATRIX_MUL_SIMD_IKJ:
			kernel_str = "mm_mul_simd_ikj";

			if (N < CACHE_CONST_BLOCKING_SIZE)
			{
				std::cerr << "Precompiler cache blocking size must be at least N" << std::endl;
				exit(1);
			}
			break;

		case MATRIX_MATRIX_MUL_OPENMP_IKJ:
			kernel_str = "mm_mul_simd_openmp_ikj";

			if (N < CACHE_CONST_BLOCKING_SIZE)
			{
				std::cerr << "Precompiler cache blocking size must be at least N" << std::endl;
				exit(1);
			}
			break;
	}
	std::cout << " + kernel: " << kernel_str << std::endl;

	/*
	 * Initialization
	 */
	switch (variant_id)
	{
		default:
			std::cerr << "***" << std::endl;
			std::cerr << "Setup not implemented" << std::endl;
			std::cerr << "***" << std::endl;
			exit(-1);
			break;

		case MATRIX_SUM_ROWWISE:
		case MATRIX_SUM_COLWISE:
			matrix_setup_A(N, A);
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IJK:
		case MATRIX_MATRIX_MUL_SIMPLE_JIK:
		case MATRIX_MATRIX_MUL_SIMPLE_IKJ:
		case MATRIX_MATRIX_MUL_SIMPLE_JKI:
		case MATRIX_MATRIX_MUL_SIMPLE_KIJ:
		case MATRIX_MATRIX_MUL_SIMPLE_KJI:
		case MATRIX_MATRIX_MUL_RESTRICTED_IKJ:
		case MATRIX_MATRIX_MUL_VAR_BLOCKED_IKJ:
		case MATRIX_MATRIX_MUL_BLOCKED_IKJ:
		case MATRIX_MATRIX_MUL_SIMD_IKJ:
		case MATRIX_MATRIX_MUL_OPENMP_IKJ:
			matrix_setup_A(N, A);
			matrix_setup_B(N, B);
			matrix_zero_C(N, C);
			break;
	}


	std::cout << std::endl;
	std::cout << "Starting benchmark... " << std::flush;

	/**
	 * Return buffer for scalar
	 */
	double retscalar;

	Stopwatch stopwatch;

	// Number of iterations
	int num_iterations = 0;
	{
		/*
		 * Run benchmark for moderate problem sizes
		 */
		if (N < 512)
			retscalar = run_benchmark(variant_id, N, A, B, C, cache_blocking_size);

		while (true)
		{
			matrix_zero_C(N, C);
			flush_cache();

			/**
			 * Run matrix multiplication
			 */
			stopwatch.start();

			retscalar = run_benchmark(variant_id, N, A, B, C, cache_blocking_size);

			stopwatch.stop();

			num_iterations++;

			// Stop iterations after 1 second(s)
			if (stopwatch() > 1)
				break;

			// Or after 10 iterations
			if (num_iterations >= 10)
				break;
		}
	}

	std::cout << "Finished" << std::endl;
	std::cout << std::endl;

	double elapsed_time = stopwatch();

	double num_flops = -1;


	if (N <= 8)
	{
		std::cout << "Matrix A:" << std::endl;
		print_matrix(N, A);
		std::cout << std::endl;

		std::cout << "Matrix B:" << std::endl;
		print_matrix(N, B);
		std::cout << std::endl;

		std::cout << "Matrix C:" << std::endl;
		print_matrix(N, C);
		std::cout << std::endl;
	}


	/**
	 * Postprocessing (validation, etc.)
	 */
	switch (variant_id)
	{
		default:
			std::cerr << "Validation not implemented" << std::endl;
			exit(-1);
			break;

		case MATRIX_SUM_ROWWISE:
		case MATRIX_SUM_COLWISE:
			{
				double retscalar_test = 0;
				for (std::size_t i = 0; i < N; i++)
					for (std::size_t j = 0; j < N; j++)
						retscalar_test += A[i*N+j];

				if (std::abs(retscalar-retscalar_test) > 1e-8)
				{
					std::cout << "retvalue: " << retscalar << std::endl;
					std::cout << "should be: " << retscalar_test << std::endl;
					std::cerr << "Sanity check failed!" << std::endl;
					exit(-1);
				}
			}
			// each matrix element is added to the scalar
			num_flops = N*N;
			break;

		case MATRIX_MATRIX_MUL_SIMPLE_IJK:
		case MATRIX_MATRIX_MUL_SIMPLE_JIK:
		case MATRIX_MATRIX_MUL_SIMPLE_IKJ:
		case MATRIX_MATRIX_MUL_SIMPLE_JKI:
		case MATRIX_MATRIX_MUL_SIMPLE_KIJ:
		case MATRIX_MATRIX_MUL_SIMPLE_KJI:

		case MATRIX_MATRIX_MUL_RESTRICTED_IKJ:
		case MATRIX_MATRIX_MUL_VAR_BLOCKED_IKJ:
		case MATRIX_MATRIX_MUL_BLOCKED_IKJ:

		case MATRIX_MATRIX_MUL_SIMD_IKJ:
		case MATRIX_MATRIX_MUL_OPENMP_IKJ:
			validate_matrix_C(N, C);

			// for each output element of C (N*N elements), there are (N multiplications and N-1 additions)
			num_flops = N*N*(N+N-1);
			break;
	}


	/**
	 * Output information
	 */
	std::cout << " + elapsed_time: " << elapsed_time << std::endl;
	std::cout << " + num_iterations: " << num_iterations << std::endl;
	std::cout << " + time/iteration: " << elapsed_time/num_iterations << std::endl;
	std::cout << " + num_flops: " << num_flops << std::endl;
	std::cout << " ++ g_num_flops: " << num_flops*1e-9 << std::endl;
	std::cout << " ++ t_num_flops: " << num_flops*1e-12 << std::endl;
	std::cout << " ++ g_num_flops/s: " << num_flops*1e-9/elapsed_time*num_iterations << std::endl;
	std::cout << " ++ t_num_flops/s: " << num_flops*1e-12/elapsed_time*num_iterations << std::endl;

	/**
	 * Free allocated data
	 */
	delete [] A;
	delete [] B;
	delete [] C;
}
