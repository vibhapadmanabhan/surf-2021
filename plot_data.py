import matplotlib.pyplot as plt

with open("fO2.txt", "r") as f:
    y = f.readlines()

with open("planet_size.txt", "r") as f:
    x = f.readlines()

plt.plot(x, y)
plt.show()