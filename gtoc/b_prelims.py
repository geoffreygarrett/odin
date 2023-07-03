import numpy as np
import matplotlib.pyplot as plt

# Parameters
rho = 0.004
N = np.linspace(2, 100, 500)  # Vary N from 0 to 100

# Solve the function for M (assuming rho=1)
M = np.log(N / 2) / rho

# Create the plot
plt.figure()
plt.plot(N, M)  # Swap N and M
plt.xlabel("N")
plt.ylabel(r"$\overline{M}$")
plt.title("Average mass allowed for N=2 to 100")
plt.xlim([0, 100])  # Ensure that the x-axis ranges from 0 to 100
plt.grid(True)
plt.show()
