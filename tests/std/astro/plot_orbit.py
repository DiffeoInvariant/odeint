from sys import argv
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':

    if len(argv) > 1:
        fln = argv[1]
    else:
        fln = 'orbit1.csv'

    data = pd.read_csv(fln)
    t = list(data['t'])
    x = list(data['x_0'])
    y = list(data['x_1'])
 
    print(f"Trajectory contains {len(t)} timesteps.\n")
    plt.figure()
    plt.scatter(x, y, c='black',s=0.01)
    plt.title('Integrated orbit')
    plt.show()
