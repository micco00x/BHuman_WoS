import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

# Read file composed of list of support foot at mpc_iter = 0:
def read_support_foot_file(path):
    l_supp = []
    with open(path) as f:
        for l in f.readlines():
            p = l.rstrip('\n').split()
            x = float(p[0])
            y = float(p[1])
            z = float(p[2])
            theta = float(p[3])
            l_supp.append((x, y, z, theta))
    return l_supp

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

    l_supp = read_support_foot_file(supp_file_path)

    fig = plt.figure()

    ### Subplot 1 ###
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    ax1.plot(com_trajectory.x, com_trajectory.y, label='CoM')
    ax1.plot(zmp_trajectory.x, zmp_trajectory.y, label='ZMP')
    #ax1.plot(lsole_trajectory.x, lsole_trajectory.y, label='lsole')
    #ax1.plot(rsole_trajectory.x, rsole_trajectory.y, label='rsole')

    support_foot_size = 0.05
    for support_foot_configuration in l_supp:
        (x, y, z, theta) = support_foot_configuration
        x0 = x - support_foot_size / 2.0
        y0 = y - support_foot_size / 2.0
        rect = patches.Rectangle(
            (x0, y0), support_foot_size, support_foot_size,
            linewidth=1, edgecolor='g', facecolor='none')
        ax1.add_patch(rect)

    ax1.legend()
    ax1.grid()
    ax1.axis('equal')

    ### Subplot 2 ###
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.set_xlabel('t')
    ax2.set_ylabel('z')

    delta_t = 0.01
    samples = lsole_trajectory.z.shape[0]

    tt = np.linspace(0.0, delta_t * samples, samples)
    ax2.plot(tt, lsole_trajectory.x, label='lsole x')
    ax2.plot(tt, lsole_trajectory.y, label='lsole y')
    ax2.plot(tt, rsole_trajectory.x, label='rsole x')
    ax2.plot(tt, rsole_trajectory.y, label='rsole y')
    ax2.plot(tt, lsole_trajectory.z, label='lsole')
    ax2.plot(tt, rsole_trajectory.z, label='rsole')

    ax2.legend()
    ax2.grid()

    ### Subplot 3 ###
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_xlabel('t')
    ax3.set_ylabel('ZMP.x')

    ax3.plot(tt, zmp_trajectory.x, label='ZMP.x')

    ax3.legend()
    ax3.grid()
    
    ### Subplot 4 ###
    ax3 = fig.add_subplot(2, 2, 4)
    ax3.set_xlabel('t')
    ax3.set_ylabel('ZMP.y')

    ax3.plot(tt, zmp_trajectory.y, label='ZMP.y')

    ax3.legend()
    ax3.grid()

    plt.show()