import matplotlib.pyplot as plt

# Define function to read data from files
def read_data(file_path):
    with open(file_path, 'r') as file:
        data = file.readlines()
    processors = []
    times = []
    for line in data:
        if ':' in line:
            proc, time = line.strip().split(':')
            processors.append(int(proc))
            times.append(float(time))
    return processors, times

# File paths and corresponding labels
files = {
        "nos3_2_output.txt": "nos3",
    "bcsstk27_3_output.txt": "bcsstk27",
    "bodyy5_4_output.txt": "bodyy5",
        "Kuu_5_output.txt": "Kuu",
    "Dubcova2_6_output.txt": "Dubcova2"


}

# Initialize variables to track max observed SI and best performance
max_si_observed = 0
best_wall_clock_times = {}
best_si_processors = {}

# Process each dataset to calculate max SI and best performance
si_data = {}
for file_name, label in files.items():
    file_path = f"./{file_name}"  # Adjust file path as needed
    processors, times = read_data(file_path)
    time_at_1_proc = times[0]
    si = [time_at_1_proc / time for time in times]
    si_data[label] = (processors, si)
    max_si_observed = max(max_si_observed, max(si))
    
    # Find the processor count with minimum wall-clock time
    min_time_idx = times.index(min(times))
    best_wall_clock_times[label] = (processors[min_time_idx], times[min_time_idx])
    
    # Find the processor count with maximum Speed Improvement
    max_si_idx = si.index(max(si))
    best_si_processors[label] = (processors[max_si_idx], si[max_si_idx])

# Initialize plots
plt.figure(figsize=(12, 10))

# Wall-clock time plot
plt.subplot(2, 1, 1)
for file_name, label in files.items():
    file_path = f"./{file_name}"  # Adjust file path as needed
    processors, times = read_data(file_path)
    plt.plot(processors, times, marker='.', label=label)
    
   #  # Mark the best performance (minimum wall-clock time)
   #  best_proc, best_time = best_wall_clock_times[label]
   #  plt.scatter(best_proc, best_time, color='black', zorder=5)
   #  plt.text(
   #      best_proc, best_time,
   #      f"  {best_proc}P\n  {best_time:.2f}s",
   #      color="black",
   #      fontsize=9
   #  )

plt.title("Wall-Clock Time vs Total Processors", fontsize=14)
plt.xlabel("Total Processors", fontsize=12)
plt.ylabel("Wall-Clock Time (s)", fontsize=12)
plt.legend(title="Matrices")
plt.grid(True)

# Speed Improvement plot
plt.subplot(2, 1, 2)
for label, (processors, si) in si_data.items():
    plt.plot(processors, si, marker='.', label=label)
    
   #  # Mark the processor count achieving maximum SI
   #  best_proc, best_si = best_si_processors[label]
   #  plt.scatter(best_proc, best_si, color='black', zorder=5)
   #  plt.text(
   #      best_proc, best_si,
   #      f"  {best_proc}P\n  SI={best_si:.2f}",
   #      color="black",
   #      fontsize=9
   #  )

# Add ideal speedup line
max_processors = max(max(p for p, _ in si_data.values()))
plt.plot(
    range(1, max_processors + 1),
    range(1, max_processors + 1),
    linestyle='--',
    color='red',
    label='Ideal Speed Improvement (ISI)'
)

# Scale y-axis based on observed max SI
plt.ylim(0, max_si_observed * 1.1)  # Add a 10% margin above max observed SI

# Customize Speed Improvement plot
plt.title("Speed Improvement (SI) vs Total Processors", fontsize=14)
plt.xlabel("Total Processors", fontsize=12)
plt.ylabel("Speed Improvement (SI)", fontsize=12)
plt.legend(title="Matrices")
plt.grid(True)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
