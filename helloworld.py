import numpy as np
import matplotlib.pyplot as plt

t_burn1 = 168 # first stage burn time, seconds
m_e = 5.97219e+24 # earth mass, kilograms
m_d = 496200 + 123000 + 137000 # dry mass, kilograms
m_p = 2214000 - 137000 # propellant mass, kilograms
m_w = 2214000 + 496200 + 123000 # wet mass, kilograms
v_e = 2600 # exhaust velocity, meters per second
delta_v = v_e * np.log(m_w / m_d) # total change in the rockets velocity
d = 6371000 # starting height of the rocket relative to the earth, meters
g = -9.81 # gravitational constant for earth
G = 6.6743e-11 # Gravitational constant

t_step_size = 0.0001

Q = delta_v / (v_e * t_burn1)

def mass(time):
    return m_p * np.exp(-Q * time) + m_d

def acceleration(height, time):
    left = (m_p * delta_v * np.exp(-delta_v / v_e))/ (t_burn1 * mass(time))
    right = (G * m_e) / (height**2)
    return np.abs(right - left)

def euler_method(h0, v0, a0, t0, tf, dt):

    times = np.arange(0, tf, dt)
    heights = np.zeros(len(times))
    velocities = np.zeros(len(times))
    acceleration_values = np.zeros(len(times))
    
    heights[0] = h0
    velocities[0] = v0
    acceleration_values[0] = a0
    
    for i in range(1, len(times)):
        t = times[i-1]
        
        current_acceleration = acceleration(heights[i-1], t)
        acceleration_values[i] = current_acceleration

        heights[i] = heights[i-1] + velocities[i-1] * t_step_size
        velocities[i] = velocities[i-1] + current_acceleration * t_step_size
    return times, heights, velocities, acceleration_values

times, heights, velocities, accelerations = euler_method(d, 0, 0, 0, t_burn1, t_step_size)
print(heights)
print(np.abs((heights[-1] - heights[0]) / 1609.0))
print(accelerations[0])
print(accelerations[1])

plt.figure(figsize=(10,6))

plt.subplot(2, 1, 1)
plt.plot(times, accelerations, label="Height (m)")
plt.xlabel('Time (s)')
plt.ylabel('Height (m)')
plt.grid(True)

plt.tight_layout()
plt.show()

