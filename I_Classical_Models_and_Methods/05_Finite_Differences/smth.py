import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def planck(T, lambdas):
    h = scipy.constants.Planck
    c = scipy.constants.speed_of_light
    kb = scipy.constants.k

    b = h * c / (lambdas * kb * T)

    return 8 * np.pi * h / (
        lambdas**3 * np.expm1(b)
    )

lambdas = np.arange(300, 8000, 3) * 1e-9

fig, ax = plt.subplots()

line, = ax.plot([], [])

ax.set_xlim(lambdas.min(), lambdas.max())

maximum = planck(700 + 99 * 25, lambdas).max()
ax.set_ylim(0, maximum * 1.5)

ax.set_xlabel("Wavelength (m)")
ax.set_ylabel("Intensity")

def update(frame):
    temperature = 700 + frame * 50
    intensity = planck(temperature, lambdas)

    line.set_data(lambdas, intensity)
    ax.set_title(f"T = {temperature} K")

    return line,

animation = FuncAnimation(
    fig,
    update,
    frames=100,
    interval=25,
    blit=True
)

animation.save(
    "/home/driftwave/planck_animation.gif",
    writer=PillowWriter(fps=20)
)

print("Saved to /home/driftwave/planck_animation.gif")

plt.show()