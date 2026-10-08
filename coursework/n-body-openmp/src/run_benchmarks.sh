#!/bin/bash

# Array of thread counts to test
threads=(1 2 4 8 16 24 32 48) # Adjust this if needed

# Loop through different thread counts
for num_threads in "${threads[@]}"; do
    echo "Running with $num_threads threads..."
    export OMP_NUM_THREADS=$num_threads
    start=$(date +%s%3N)
    ./main
    end=$(date +%s%3N)
    elapsed_time=$((end - start))

    echo "Time with $num_threads threads: $((elapsed_time)) ms"
    echo "$num_threads $elapsed_time" >> times.dat # Save to file
done

echo "Execution times saved to times.dat"