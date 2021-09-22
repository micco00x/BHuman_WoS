import matplotlib.pyplot as plt
import numpy as np

class Trajectory:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

# Read file composed of one position per line:
def read_position_file(path):
    # Init lists:
    xx = []
    yy = []
    zz = []

    # Read file:
    with open(path) as f:
        for l in f.readlines():
            p = l.rstrip('\n').split()
            xx.append(float(p[0]))
            yy.append(float(p[1]))
            zz.append(float(p[2]))

    # Setup np arrays:
    xx = np.array(xx)
    yy = np.array(yy)
    zz = np.array(zz)
    trajectory = Trajectory(xx, yy, zz)
    return trajectory

if __name__ == '__main__':
    com_file_path = '/tmp/com.txt'
    zmp_file_path = '/tmp/zmp.txt'
    lsole_file_path = '/tmp/lsole.txt'
    rsole_file_path = '/tmp/rsole.txt'
    supp_file_path = '/tmp/supp.txt'

    com_trajectory = read_position_file(com_file_path)
    zmp_trajectory = read_position_file(zmp_file_path)
    lsole_trajectory = read_position_file(lsole_file_path)
    rsole_trajectory = read_position_file(rsole_file_path)

    plt.subplot(2, 1, 1)
    plt.title('CoM/ZMP')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.plot(com_trajectory.x, com_trajectory.y, label='CoM')
    plt.plot(zmp_trajectory.x, zmp_trajectory.y, label='ZMP')

    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.title('Feet')
    plt.xlabel('t')
    plt.ylabel('z')

    delta_t = 0.01
    samples = lsole_trajectory.z.shape[0]

    tt = np.linspace(0.0, delta_t * samples, samples)
    plt.plot(tt, lsole_trajectory.z, label='lsole')
    plt.plot(tt, rsole_trajectory.z, label='rsole')

    plt.legend()
    plt.grid()
    plt.show()