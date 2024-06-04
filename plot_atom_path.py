import matplotlib.pyplot as plt
import numpy as np


position = np.loadtxt("atom_path.txt")

#print pos
x = position[:, 0]
y = position[:, 1]

plt.figure(figsize=(8,8))

#plt.plot(x, y, '-', linewidth=0.5)
plt.scatter(x, y, s=30 , c=pos[:,3])     #s=pos[:,3]/1000  waht is pos[:,2] pos[0,3]?...?
plt.title("Atom diffusion path")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(-11, 10)
plt.ylim(-11, 10)
plt.grid(True)
plt.show()