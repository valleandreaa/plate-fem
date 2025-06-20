import pandas as pd
import matplotlib.pyplot as plt
def PlotConvergenceOpenFOAM(pathDat='C:/Users/Admin/OneDrive/Desktop/laminar/Scenario_1_laminar/postProcessing/PressureInlet/0/surfaceFieldValue.dat', xlabel, ylabel, title ):
    df = pd.read_csv(pathDat, delim_whitespace=True,  header=None, skiprows=5)

    plt.plot(df.loc[:,0], df.loc[:,1])
    plt.title("Sports Watch Data")
    plt.xlabel("Average Pulse")
    plt.ylabel("Calorie Burnage")
    plt.show()
