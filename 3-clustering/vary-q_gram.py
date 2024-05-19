import os
import subprocess
import time
import matplotlib.pyplot as plt

# Define the different sizes for gram_size
gram_sizes = range(1, 7)  # Example: from 1 to 6

config_file = '/home/s/shaoqi/DNAStorageToolkit/3-clustering/clustering_config.cfg'
output_dir = '/home/s/shaoqi/DNAStorageToolkit/3-clustering/output'

# Function to update gram_size in the source code
def update_gram_size(size):
    file_path = '/home/s/shaoqi/DNAStorageToolkit/3-clustering/src/parallel_clustering.cpp'
    with open(file_path, 'r') as file:
        data = file.readlines()
    for i, line in enumerate(data):
        if line.startswith('const int gram_size ='):
            data[i] = f'const int gram_size = {size};\n'
            break
    with open(file_path, 'w') as file:
        file.writelines(data)

# Function to compile the source code
def compile_code():
    result = subprocess.run(['make', '-C', '/home/s/shaoqi/DNAStorageToolkit/3-clustering', 'exe/clustering'])
    if result.returncode != 0:
        raise Exception('Compilation failed')

# Function to parse the output and collect metrics
def parse_output(output):
    lines = output.split('\n')
    metrics = {}
    for line in lines:
        if 'Accuracy (gamma=' in line:
            metrics['accuracy'] = float(line.split('Accuracy (gamma=')[1].split(') :')[1])
        elif '--- ' in line and 'seconds for' in line:
            time_taken = float(line.split('--- ')[1].split(' seconds')[0])
            if 'complete program execution' in line:
                metrics['total_time'] = time_taken
    return metrics

# Function to run clustering and collect metrics
def run_clustering_and_collect_metrics(size):
    update_gram_size(size)
    compile_code()

    start_time = time.time()
    
    # Run clustering process
    result = subprocess.run(['make', '-C', '/home/s/shaoqi/DNAStorageToolkit/3-clustering', 'run'], capture_output=True, text=True)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # Collect metrics from the output
    metrics = parse_output(result.stdout)
    metrics['gram_size'] = size
    metrics['total_time'] = elapsed_time
    return metrics

# Function to print metrics in a formatted manner
def print_metrics(metrics):
    print("\n" + "="*50)
    print(f"Gram Size: {metrics['gram_size']}")
    print(f"Accuracy: {metrics.get('accuracy', 'N/A')}")
    print(f"Total Time: {metrics.get('total_time', 'N/A')} seconds")
    print("="*50 + "\n")

# List to store performance metrics
performance_metrics = []

# Run clustering for each gram_size and collect metrics
for size in gram_sizes:
    metrics = run_clustering_and_collect_metrics(size)
    performance_metrics.append(metrics)
    print_metrics(metrics)

# Analyze and plot the results
gram_sizes = [metrics['gram_size'] for metrics in performance_metrics]
accuracies = [metrics.get('accuracy', 0) for metrics in performance_metrics]
total_times = [metrics.get('total_time', 0) for metrics in performance_metrics]

plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(gram_sizes, accuracies, marker='o')
plt.title('Effect of Gram Size on Clustering Performance')
plt.ylabel('Accuracy')

plt.subplot(2, 1, 2)
plt.plot(gram_sizes, total_times, marker='o')
plt.ylabel('Total Time (s)')
plt.xlabel('Gram Size')

plt.tight_layout()
plt.savefig('clustering_performance_analysis.png')
plt.show()
