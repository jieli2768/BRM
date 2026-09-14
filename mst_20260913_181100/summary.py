import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import pandas as pd

# 1. Read the CSV file
df = pd.read_csv('cat.csv')

# 2. Find the largest slot value for each student
max_slots = df.groupby('student')['slot'].max()

# 3. Count how many students reached each maximum slot
slot_counts = max_slots.value_counts().sort_index()
print(f'Mean: {max_slots.mean():.2f}')
print(f'Std Dev: {max_slots.std():.2f}')

# 4. Plot the histogram / count chart
plt.figure(figsize=(8, 5))
plt.bar(
    slot_counts.index,
    slot_counts.values,
    color='skyblue',
    edgecolor='black',
    align='center',
)
plt.xlabel('Test Length')
plt.ylabel('Count of Students')
plt.title('Distribution of Test Length per Student')
plt.xticks(slot_counts.index)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

plt.show()