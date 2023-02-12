# An implementation of an Ising spin-glass of size NxN
# With fixed boundary conditions using Metropolis-Hastings
# Connectivity is initialized as a Gaussian distribution N(0, s^2/N)
# Updates occur at randomly selected sites

import numpy as np
import matplotlib.pyplot as plt

# Fix random seed
np.random.seed(123)

# Set size of model N and initial spins
N = 32
spins = np.random.choice([-1, 1], (N, N))

# Fix number of timesteps and some containers
timesteps = 10000
mag = np.zeros(timesteps+1)
energy = np.zeros(timesteps+1)

# Initialize interaction arrays
# Have 4 arrays: up, down, left right
# These represent the interaction strengths to the
# up/down/left/right neighbours of a site
# There is a symmetry between these matrices
# This is not the most memory efficient solution
s_h = 1
s_v = 1
up = np.zeros((N, N))
down = np.zeros((N, N))
left = np.zeros((N, N))
right = np.zeros((N, N))

up[1:N, :] = np.random.rand(N-1, N) * s_v
down[0:N-1, :] = up[1:N, :]
left[:, 1:N] = np.random.rand(N, N-1) * s_h
right[:, 0:N-1] = left[:, 1:N]

mag[0] = spins.sum()

for i in range(N):
    for j in range(N):
        if i == 0:
            up_neighbour = 0
            down_neighbour = spins[i+1, j]
        elif i == N-1:
            up_neighbour = spins[i-1, j]
            down_neighbour = 0
        else:
            up_neighbour = spins[i-1, j]
            down_neighbour = spins[i+1, j]
        if j == 0:
            left_neighbour = 0
            right_neighbour = spins[i, j+1]
        elif j == N-1:
            left_neighbour = spins[i, j-1]
            right_neighbour = 0
        else:
            left_neighbour = spins[i, j-1]
            right_neighbour = spins[i, j+1]

        energy[0] += spins[i, j]*(up[i, j]*up_neighbour + down[i, j]*down_neighbour +
                                  left[i, j]*left_neighbour + right[i, j]*right_neighbour)

# Avoid double count - each neighbour pair
# counted twice in above since loop over each site
energy[0] /= 2

# Fix beta (inverse temerature) - from analysis we know that
# system in glassy-phase for T<s so beta>1/s. Performance
# of random updates isn't good so don't select temperature
# too low
beta = 1/4

# Define proposal step


def proposal(s_array):
    _N = s_array.shape[0]
    return np.random.choice(_N, 2)


def energy_change(spin_site, bt, s_array, up_array, down_array, left_array, right_array):
    i = spin_site[0]
    j = spin_site[1]

    if i == 0:
        up_neighbour = 0
        down_neighbour = s_array[i+1, j]
    elif i == N-1:
        up_neighbour = s_array[i-1, j]
        down_neighbour = 0
    else:
        up_neighbour = s_array[i-1, j]
        down_neighbour = s_array[i+1, j]
    if j == 0:
        left_neighbour = 0
        right_neighbour = s_array[i, j+1]
    elif j == N-1:
        left_neighbour = s_array[i, j-1]
        right_neighbour = 0
    else:
        left_neighbour = s_array[i, j-1]
        right_neighbour = s_array[i, j+1]

    dE_tmp = 2*s_array[i, j]*(up_array[i, j]*up_neighbour + down_array[i, j]*down_neighbour +
                              left_array[i, j]*left_neighbour + right_array[i, j]*right_neighbour)
    return dE_tmp


def acceptance(bt, energy):
    if energy <= 0:
        return -1
    else:
        prob = np.exp(-bt*energy)
        if prob > np.random.random():
            return -1
        else:
            return 1


# Define update step
dE = 0
dM = 0


def update(bt, s_array, up_array, down_array, left_array, right_array):
    global dE
    global dM

    # Proposal Step
    site = proposal(s_array)

    # Calculate energy change
    dE = energy_change(site, bt, s_array, up_array,
                       down_array, left_array, right_array)
    dM = -2*s_array[site[0], site[1]]

    # Acceptance step
    accept = acceptance(bt, dE)

    if accept == -1:
        s_array[site[0], site[1]] *= -1
    else:
        dE = 0
        dM = 0

    return s_array


def _main_loop(ts, s_array, up_array, down_array, left_array, right_array):
    s_temp = s_array.copy()
    for i in range(ts):
        update_step = update(beta, s_temp, up_array,
                             down_array, left_array, right_array)
        s_temp = update_step
        energy[i+1] = energy[i] + dE
        mag[i+1] = mag[i] + dM



def _main_loop_SA(ts, bt_initial, s_array, up_array, down_array, left_array, right_array):
    s_temp = s_array.copy()
    bt_live = bt_initial
    for i in range(ts):
        if ts % 500 == 0:
            bt_live *= 1/0.9 #cooling schedule
        update_step = update(bt_live, s_temp, up_array,
                             down_array, left_array, right_array)
        s_temp = update_step
        energy[i+1] = energy[i] + dE
        mag[i+1] = mag[i] + dM


#### Run Main Loop
_main_loop_SA(timesteps, beta, spins, up, down, left, right)
mag = mag / (N**2)
energy = energy / (N**2)

# plot magnetism and energy evolving in time
fig, ax1 = plt.subplots()
ax1.set_xlabel("Time step")
ax1.set_ylabel("Magnetism", color='blue')
ax1.plot(mag, color='blue')

ax2 = ax1.twinx()
ax2.set_ylabel("Energy", color='red')
ax2.plot(energy, color='red')

plt.show()
