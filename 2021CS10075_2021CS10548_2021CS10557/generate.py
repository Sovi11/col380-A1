import random
from sys import argv
# I want to take the input from the user argument
n = int(argv[1])
for i in range(n):
    for j in range(n):
        print(random.randint(-1000, 1000) , end =" ")
    print("")
