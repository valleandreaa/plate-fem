import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

def PlotError(x, y, Xlabel, Ylabel, NameFile):
    plt.figure(1)

    plt.yscale('log')
    plt.xscale('log')
    plt.plot(x, y)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.tight_layout()
    # plt.legend()
    plt.savefig(NameFile)
x=[35, 234, 985, 3941]
y=[135.726, 129.073, 124.149, 121.975]

PlotError(x, y, Xlabel='DOF', Ylabel='Norm of the error in the stress', NameFile='./error_dof.png')
