import matplotlib.pyplot as plt

def PlotDisplacement(x, y, i, Xlabel, Ylabel, NameFile):
    plt.figure(1)
    plt.plot(x, y, label=i)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.legend()
    plt.savefig(NameFile)

def PlotError(x, y, Xlabel, Ylabel, NameFile):
    plt.figure(2)

    plt.yscale('log')
    plt.xscale('log')
    plt.plot(x, y)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.legend()
    plt.savefig(NameFile)

def PlotErrorAlpha(x, y, label, Xlabel, Ylabel, NameFile):
    plt.figure(3)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(x, y, label=label)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.legend()
    plt.savefig(NameFile)