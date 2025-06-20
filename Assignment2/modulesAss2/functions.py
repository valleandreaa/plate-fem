import math
from numpy.linalg import multi_dot
import numpy as np

def AnalyticalFuncions(mesh,nu, E, sigma,a,  r, tetha):
    """
    Analytical Functions
    mesh: instance mesh
    nu: Poisson's ratio
    E: Young's modulus
    sigma: Applied tension
    a: parameter "a" of an hole in the plate
    r: radius of an hole in the plate
    tetha: angle of the point
    """

    nu_=nu/(1-nu)
    E_=E/(1-nu**2)

    Ux=((1+nu_)*sigma/E_)*((1/(1+nu_))*r*math.cos(tetha)+
                          (2 / (1 + nu_))*a**2/r*math.cos(tetha)+
                          0.5*a**2/r*math.cos(3*tetha)-
                          0.5*a**4/r**3*math.cos(3*tetha))

    Uy=((1+nu_)*sigma/E_)*((-nu_/(1+nu_))*r*math.sin(tetha)-
                          ((1-nu_) / (1 + nu_))*a**2/r*math.sin(tetha)+
                          0.5*a**4/r**3*math.cos(3*tetha))

    Sigma_xx=sigma*(1-a**2/r**2*(1.5*math.cos(2*tetha)+math.cos(4*tetha))+3*a**4*math.cos(4*tetha)/(2*r**4))
    Sigma_yy=sigma*(-a**2/r**2*(0.5*math.cos(2*tetha)-math.cos(4*tetha))-3*a**4*math.cos(4*tetha)/(2*r**4))
    Sigma_xy=sigma*(-a**2/r**2*(0.5*math.sin(2*tetha)+math.sin(4*tetha))+3*a**4*math.sin(4*tetha)/(2*r**4))

    return Ux, Uy, Sigma_xx, Sigma_yy, Sigma_xy

def EnforcingSymDispl(mesh,xCoord, yCoord):
    E = 2.1 * pow(10, 7)
    nu = 0.29
    sigma = 100
    a = 2.5


    r = math.sqrt(xCoord ** 2 + yCoord ** 2)
    tetha = math.atan(yCoord / xCoord)

    Ux, Uy, Sigma_xx, Sigma_yy, Sigma_xy = AnalyticalFuncions(mesh, nu, E, sigma, a, r, tetha)
    return Ux,Uy
def ExactCalc(mesh,e):
    """
    Calculation of the point of the cell
    mesh: instance of mesh
    e: number of the element
    """
    E = 2.1*pow(10,7)
    nu = 0.29
    sigma=100
    a=2.5

    xCoord=mesh.points[mesh.elements[e][0]][0]*1/3+mesh.points[mesh.elements[e][1]][0]*1/3+mesh.points[mesh.elements[e][2]][0]*1/3
    yCoord=mesh.points[mesh.elements[e][0]][1]*1/3+mesh.points[mesh.elements[e][1]][1]*1/3+mesh.points[mesh.elements[e][2]][1]*1/3

    r =math.sqrt(xCoord**2+yCoord**2)
    tetha=math.atan(yCoord/xCoord)

    Ux, Uy, Sigma_xx, Sigma_yy, Sigma_xy = AnalyticalFuncions(mesh, nu, E, sigma, a, r, tetha)
    VonMisesExact=math.sqrt(Sigma_xx ** 2 + Sigma_yy ** 2 - Sigma_xx * Sigma_yy +  3 * Sigma_xy**2)


    return Ux, Uy, Sigma_xx, Sigma_yy, Sigma_xy, VonMisesExact
def error(error, mesh):
    """
    Norm of the error in the stress
    error: error of Von Mises
    mesh: instance of mesh
    """
    error = np.array(error)
    err=0
    for e in range(len(mesh.elements)):
        pointsEl=mesh.points
        element_coordinates = mesh.elements[e]
        n1 = element_coordinates[0]
        n2 = element_coordinates[1]
        n3 = element_coordinates[2]
        x1 = pointsEl[n1, 0]
        y1 = pointsEl[n1, 1]
        x2 = pointsEl[n2, 0]
        y2 = pointsEl[n2, 1]
        x3 = pointsEl[n3, 0]
        y3 = pointsEl[n3, 1]
        y23 = y2 - y3
        y31 = y3 - y1
        y12 = y1 - y2
        x32 = x3 - x2
        x13 = x1 - x3
        x21 = x2 - x1
        x32 = x3 - x2
        y13 = y1 - y3
        x23 = x2 - x3
        detJ = x13 * y23 - y13 * x23

        err += (0.5*detJ*(error[e])*(error[e]))
    err= err**0.5
    return err

