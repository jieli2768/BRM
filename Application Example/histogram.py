import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14   # axis label size
plt.rcParams['xtick.labelsize'] = 12  # X-axis tick numbers
plt.rcParams['ytick.labelsize'] = 12  # Y-axis tick numbers
plt.rcParams['axes.titlesize'] = 20  # Subplot titles ('slot 0', etc.)
plt.rcParams['figure.titlesize'] = 20  # Main figure title ('Bayesian Estimation')


df = pd.read_csv('cat.csv')

for student in [1,2]:
    
    df_student = df[df['student'] == student]

    # 2. Extract v1 through v5 based on their respective slots
    v1 = df_student.loc[df_student['slot'] == 0, 'bayesian_draws'].values[0]
    v2 = df_student.loc[df_student['slot'] == 10, 'bayesian_draws'].values[0]
    v3 = df_student.loc[df_student['slot'] == 20, 'bayesian_draws'].values[0]
    v4 = df_student.loc[df_student['slot'] == 30, 'bayesian_draws'].values[0]
    v5 = df_student.loc[df_student['slot'] == 40, 'bayesian_draws'].values[0]


    truetheta = df_student.loc[df_student['slot'] == 0, 'true_theta'].values[0]

    vars_dict = {
        'Position 0': v1, 
        'Position 10': v2, 
        'Position 20': v3, 
        'Position 30': v4, 
        'Position 40': v5
    }

    # 3. Parse all data first to find the global minimum and maximum for a unified x-axis range
    parsed_data = {}
    all_values = []

    for name, val_str in vars_dict.items():
        data = np.array([float(x) for x in val_str.split(';') if x.strip() != ''])
        parsed_data[name] = data
        all_values.extend(data)

    global_min = -3
    global_max = 3

    # Add a tiny padding (5%) so the bars don't touch the outer borders
    padding = (global_max - global_min) * 0.05
    if padding == 0:
        padding = 1  
    x_limits = (global_min - padding, global_max + padding)

    # Define consistent bin edges across the entire global range so bar widths match perfectly
    shared_bins = np.linspace(x_limits[0], x_limits[1], 61)

    # 4. Create a 1x5 layout figure with shared y-axis and unified x-limits
    fig, axes = plt.subplots(1, 5, figsize=(25, 5), sharey=True)

    mytitle = 'Student B (Low Ability)' if student==1 else 'Student A (High Ability)'
    fig.suptitle(mytitle, fontweight='bold', y=1.05)

    for ax, (name, data) in zip(axes, parsed_data.items()):
        # Calculate statistics
        mean_val = np.mean(data)
        sd_val = np.std(data, ddof=1) # Sample standard deviation
        
        # Draw count histogram with consistent bins and styling
        ax.hist(data, bins=shared_bins, color='steelblue', edgecolor='white')
        
        # Add vertical line for truetheta on every subplot
        ax.axvline(x=truetheta, color='crimson', linestyle='--', linewidth=2, label=f'True Theta ({truetheta})')
        
        # Subplot labels and settings
        ax.set_title(name, fontweight='bold')
        ax.set_xlabel('Theta')
        ax.set_xlim(x_limits)
        
        # Force y-axis tick labels to show on every subplot
        ax.tick_params(labelleft=True)
        
        # Add text box containing Mean and SD in the top-right corner
        text_str = f"Mean: {mean_val:.2f}\nSD: {sd_val:.2f}"
        props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
        ax.text(0.95, 0.95, text_str, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', horizontalalignment='right', bbox=props)

    # Add a shared legend on the first subplot (or globally) if desired
    axes[0].legend(loc='upper left', fontsize=10)

    # Adjust layout for neat spacing
    plt.tight_layout()
    plt.show()

