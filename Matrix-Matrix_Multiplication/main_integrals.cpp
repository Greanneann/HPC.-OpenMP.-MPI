#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <cmath>

double pi_approx(const double& x){
    return 4 / (1 + x*x);
}

double dilogarithm_function(const double& x){
     return - std::log(1 - x) / x;
}

double hand_written_reduction(const int& num_threads, const int& num_intervals, double (*func)(const double&)) {

    double start_time, end_time;

    // Set precision for output
    std::cout << std::setprecision(8);
    std::cerr << std::setprecision(8);

    // Record the start time
    start_time = omp_get_wtime();

    // Set the number of threads
    omp_set_num_threads(num_threads);

    // Output thread count and number of intervals if verbose is enabled
    std::cout << "---------------------------------------\n";
    std::cout << " + num_threads: " << num_threads << std::endl;
    std::cout << " + num_intervals: " << num_intervals << std::endl;
    std::cout << "---------------------------------------\n";

    double weight = 1.0 / num_intervals;
    double sum = 0.0;

    #pragma omp parallel
    {  
        // Private variable for partial sum
        double partial_sum = 0.0;

        // Parallel for loop
        #pragma omp for
        for (long i = 0; i < num_intervals; i++) {
            double x = weight * ((double)i + 0.5); // Compute x for the current interval
            double y = func(x);                   // Apply the provided function
            partial_sum += weight * y;            // Accumulate partial result
        }

        #pragma omp critical
        std::cout << "Thread " << omp_get_thread_num() << " partial sum = " << partial_sum << std::endl;
        sum += partial_sum;
    }

    // Check if the approximation is precise enough
    // if (std::abs(sum - 3.14159265358979323846) <= 1e-8){
    //     std::cout << "Close enough!\n";
    // }

    // Output final results
    std::cout << "---------------------------------------\n";
    std::cout << "Approximation: " << sum << std::endl;

    // Record the end time
    end_time = omp_get_wtime();

    // Output elapsed time

    std::cout << "---------------------------------------\n";
    printf("Elapsed time: %f seconds\n", end_time - start_time);

    // Return the elapsed time
    return end_time - start_time;
}

double omp_reduction(const int& num_threads, const int& num_intervals, double (*func)(const double&)) {
    double start_time, end_time;

    // Set precision for output
    std::cout << std::setprecision(8);
    std::cerr << std::setprecision(8);

    // Record the start time
    start_time = omp_get_wtime();

    // Set the number of threads
    omp_set_num_threads(num_threads);

    // Output thread count and number of intervals if verbose is enabled
    std::cout << "---------------------------------------\n";
    std::cout << " + num_threads: " << num_threads << std::endl;
    std::cout << " + num_intervals: " << num_intervals << std::endl;
    std::cout << "---------------------------------------\n";

    double weight = 1.0 / num_intervals;
    double sum = 0.0;

    // Parallel for loop with reduction
    #pragma omp parallel for reduction(+:sum)
    for (long i = 0; i < num_intervals; i++) {
        double x = weight * ((double)i + 0.5); // Compute x for the current interval
        double y = func(x);                   // Apply the provided function
        sum += weight * y;            // Accumulate partial result
    }
    
    // Check if the approximation is precise enough
    // if (std::abs(sum - 3.14159265358979323846) <= 1e-8){
    //     std::cout << "Close enough!\n";
    // }

    // Output final results
    std::cout << "---------------------------------------\n";
    std::cout << "Approximation: " << sum << std::endl;

    // Record the end time
    end_time = omp_get_wtime();

    // Output elapsed time

    std::cout << "---------------------------------------\n";
    printf("Elapsed time: %f seconds\n", end_time - start_time);

    // Return the elapsed time
    return end_time - start_time;
}

void benchmark(const std::string& output_file, double (*reduction)(const int&, const int&, double (*)(const double&)), double(*func)(const double&)) {

    // Define the number of threads and interval sizes to test
    std::vector<int> thread_counts = {1, 2, 4, 8, 16}; // Number of threads to test
    std::vector<long int> interval_sizes = {10000, 100000, 1000000, 10000000}; // Number of intervals for integration

    // Open the output file for writing results
    std::ofstream file(output_file);
    if (!file.is_open()) {
        // Error handling if the file cannot be opened
        std::cerr << "Error: Unable to open file " << output_file << " for writing." << std::endl;
        return;
    }

    // Write the header line with the interval sizes
    file << "num_threads/num_intervals"; // First column header
    for (auto num_intervals : interval_sizes) {
        file << "\t" << num_intervals; // Column headers for each interval size
    }
    file << std::endl;

    // Perform benchmarks for different thread counts
    for (auto num_threads : thread_counts) {
        file << num_threads; // Write the number of threads in the first column

        // Test for different interval sizes
        for (auto num_intervals : interval_sizes) {
            // Call the reduction function and measure elapsed time
            double elapsed_time = reduction(num_threads, num_intervals, func);

            // Write the elapsed time to the file
            file << "\t" << std::fixed << std::setprecision(9) << elapsed_time;
        }

        // End the current line after testing all interval sizes for the current thread count
        file << std::endl;
    }

    // Close the output file
    file.close();

    // Notify the user that the results have been saved
    std::cout << "Benchmark results saved to " << output_file << std::endl;
}

int main(int argc, char *argv[])
{
    // Use pi_approx for the calculation of pi.
    // Use the dilogarithm_function to calculate Li2(x) at x = 1.

    std::string output_file_hand_written_reduction = "hand_written_reduction.csv";
    benchmark(output_file_hand_written_reduction, hand_written_reduction, dilogarithm_function);

    // std::string output_file_omp_reduction = "omp_reduction.csv";
    // benchmark(output_file_omp_reduction, omp_reduction, dilogarithm_function);
    
    return 0;
}