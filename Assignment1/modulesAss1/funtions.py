import math
import matplotlib.pyplot as plt
def ErrorFunction(err,hSize,h, alpha, P, epsilon_h):
    '''
    Calculation of the error
    err: list of errors
    hSize: list of size of elements
    alpha: coefficient
    P: force applied
    epsilon_h: error from the problem
    '''
    epsilon= P**2/(2*alpha)


    err.append(math.sqrt(abs((epsilon-epsilon_h)/epsilon)))
    hSize.append(h)
    return err, hSize
