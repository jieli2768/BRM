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


df0 = pd.read_csv('cat.csv')

for student in [1,2]:
    
    df = df0[df0['student'] == student]
    
    fig, ax = plt.subplots(figsize=(18, 12))

    # 1. Plot the main Bayesian mean line
    ax.plot(
        df["slot"],
        df["bayesian_mean"],
        color="blue",
        linewidth=2,
        label="Posterior Mean",
    )

    t = 1  # 1 or 2
    
    # 2. Plot upper and lower bounds (optional, if you want visible border lines)
    ax.plot(
        df["slot"],
        df["bayesian_mean"] + t * df["bayesian_sd"],
        color="lightblue",
        linestyle="--",
        linewidth=1,
        label="Posterior Mean ± "+ ("1 " if t==1 else str(t)+"*")  +"SD",
    )
    ax.plot(
        df["slot"],
        df["bayesian_mean"] - t * df["bayesian_sd"],
        color="lightblue",
        linestyle="--",
        linewidth=1,
    )

    # 3. Add the shaded region between the second and third lines
    ax.fill_between(
        df["slot"],
        df["bayesian_mean"] - t * df["bayesian_sd"],
        df["bayesian_mean"] + t * df["bayesian_sd"],
        color="blue",
        alpha=0.2,  # Transparency of the shading
    )
    
    # 3. Add a horizontal line for true_theta
    # (Using .iloc[0] grabs the first value assuming it's a constant column)
    theta_val = df["true_theta"].iloc[0]
    
    ax.axhline(
        y=theta_val,
        color="red",
        linestyle="--",
        linewidth=1.5,
    )
    
    ax.text(
        x=df["slot"].iloc[0],
        y=theta_val,
        s=f" True Theta: {theta_val:.2f}",  # Format to 2 decimal places
        color="red",
        fontsize=12,
        verticalalignment="bottom",  # Positions text just above the line
        horizontalalignment="left",
    )



    # Customize labels and aesthetics
    ax.set_xlabel("Item Position")
    ax.set_ylabel("Theta")
    mytitle = 'Student B (Low Ability)' if student==1 else 'Student A (High Ability)'
    ax.set_title(mytitle)
    ax.legend(loc="upper left")
    ax.grid(True, linestyle=":", alpha=0.6)

    # Display the chart
    plt.show()

    
    
    
    