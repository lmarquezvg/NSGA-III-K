"""
Save approximation set.
"""

import numpy as np
import matplotlib.pyplot as plt

def saveApproximationSet(A, algorithm, problem, run, mode='save_all'):
    """Draws and saves a given approximation set"""
    N,m = np.shape(A)
    if (m == 2):
        plt.scatter(A[:,0],A[:,1],c='blue')
        plt.xlabel('$f_1$',rotation=0)
        plt.ylabel('$f_2$',rotation=0)
    elif (m == 3):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(30,45)
        ax.scatter(A[:,0],A[:,1],A[:,2],c='blue')
        ax.xaxis.set_rotate_label(False)
        ax.yaxis.set_rotate_label(False)
        ax.zaxis.set_rotate_label(False)
        ax.set_xlabel('$f_1$')
        ax.set_ylabel('$f_2$')
        ax.set_zlabel('$f_3$')
    else:
        for i in range(0,N):
            plt.plot(A[i])
        plt.xlabel('Objective function',rotation=0)
        plt.ylabel('Objective value')
        x = []
        labels = []
        for i in range(0,m):
            x.append(i)
            labels.append(str(i+1))
        plt.xticks(x,labels)
    plt.title(algorithm+' on '+problem)
    plt.tight_layout()
    if (mode == 'save_all'):
        np.savetxt('Results/'+algorithm+'_'+problem+'_{0:0=2d}D'.format(m)+'_R{0:0=2d}'.format(run)+'.pof',A,'%.6f',header=str(N)+' '+str(m))
        plt.savefig('Results/'+algorithm+'_'+problem+'_{0:0=2d}D'.format(m)+'_R{0:0=2d}'.format(run)+'.png')
        plt.close()
    elif (mode == 'save_fig'):
        plt.savefig('Results/'+algorithm+'_'+problem+'_{0:0=2d}D'.format(m)+'_R{0:0=2d}'.format(run)+'.png')
        plt.close()
    elif (mode == 'save_txt'):
        np.savetxt('Results/'+algorithm+'_'+problem+'_{0:0=2d}D'.format(m)+'_R{0:0=2d}'.format(run)+'.pof',A,'%.6f',header=str(N)+' '+str(m))
        plt.close()
    elif (mode == 'plot'):
        plt.show()
        plt.close()
    return
