import matplotlib.pyplot as plt
import statistics

file = open("pr76_hill.txt","r")
date = list()
for i in range(0,20):
    number = file.readline().split("\n")[0]
    date.append(float(number))

y = date
x = range(0,20)

# naming the x axis
plt.xlabel('Restarts')
# naming the y axis
plt.ylabel('Value')

# giving a title to my graph
plt.title('Pr76')

z = list()
for i in range(0,20):
    z.append(108159)

plt.plot(x ,y, label='hcb min')
plt.plot(x ,z, label='min')
plt.legend()
plt.show()