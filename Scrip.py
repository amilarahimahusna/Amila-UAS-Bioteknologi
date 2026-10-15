#open Anaconda powershell prompt (miniconda3) 
#instal numpy scripy matlplotlib 
pip install numpy scipy matplotlib

#Create and open a new Python file by typing the following command in your terminal
notepad simulasi_ugm.py

#When the dialog box appears asking to create a new file, click Yes.
#Copy and paste the complete simulation script below into your Notepad (1-7)

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 1. Kinetic Parameters
V1_max = 5.0
Km1 = 2.0
Ki = 3.0
X = 10.0
k2 = 1.0
k3 = 0.8
k4 = 0.3

# 2. Ordinary Differential Equations (ODEs) Formulation
def metabolic_system(t, y):
    A, B, P = y
    
    # Reaction velocity v1 with Allosteric Non-Competitive Feedback Inhibition by P
    v1 = (V1_max * X) / ((Km1 + X) * (1.0 + (P / Ki)))
    
    # Mass action kinetic rates (v2, v3, v4)
    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A
    
    # Mass balance differential equations
    dAdt = v1 - v2 - v4
    dBdt = v2 - v3
    dPdt = v3
    
    return [dAdt, dBdt, dPdt]

# 3. Initial Conditions and Simulation Timeframe (0 to 48 Hours)
initial_state = [0.0, 0.0, 0.0]  
t_span = (0, 48)                 
t_eval = np.linspace(0, 48, 1000) 

# 4. Numerical Integration Using RK45 Solver
solution = solve_ivp(metabolic_system, t_span, initial_state, t_eval=t_eval, method='RK45')

# 5. Figure Generation and Plot Customization
plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='[A] Metabolite A', color='royalblue', linewidth=2)
plt.plot(solution.t, solution.y[1], label='[B] Metabolite B', color='darkorange', linewidth=2)
plt.plot(solution.t, solution.y[2], label='[P] Product P', color='forestgreen', linewidth=2)

plt.title('Dynamic Simulation of Branch-Point Metabolic Pathway', fontsize=14, fontweight='bold')
plt.xlabel('Time (Hours)', fontsize=12)
plt.ylabel('Concentration (Arbitrary Units)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=11)
plt.tight_layout()

# 6. Save Visualization as a High-Resolution Image
output_image = 'metabolic_simulation_results.png'
plt.savefig(output_image, dpi=300)

# 7. Print Final Simulation Summary to Terminal
print("\n" + "="*50)
print("          SIMULATION SUCCESSFULLY COMPLETED          ")
print("="*50)
print(f"Visualization plot saved as : {output_image}")
print("-"*50)
print(f"Final Concentration Dynamics at t = 48 Hours:")
print(f" * Metabolite A : {solution.y[0][-1]:.4f}")
print(f" * Metabolit B : {solution.y[1][-1]:.4f}")
print(f" * Product P    : {solution.y[2][-1]:.4f}")
print("="*50 + "\n")

#Press Ctrl + S to save the file, then close the Notepad window.
#Go back to your Anaconda terminal and run the simulation using
python simulasi_ugm.py

#To locate and view your generated high-resolution plot
#open the local working directory using 
explorer .





