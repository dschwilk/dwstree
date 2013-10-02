import numpy

def linear_reg_fixed_intercept(x,y,y0):
    """Calculate least squares slope for regression with fixed y intercept =
    y0. x and y should be numpy arrays, y0 should be a scalar."""
    xm = numpy.vstack([x, numpy.zeros(len(x))]).T  # set up coefficient matrix for numpy.linalg.lstsq
    ynew = y-y0
    coeffs,residuals,rank,s = numpy.linalg.lstsq(xm,ynew)
#    print coeffs
    return coeffs[0]  # just return slope since that is all we estimated


def linear_reg(x, y):
    """ returns slope, intercept and r-squared as tuple"""
    coeffs = numpy.polyfit(x, y, 1)

    # r-squared
    p = numpy.poly1d(coeffs)
    # fit values, and mean
    yhat = [p(z) for z in x]
    ybar = sum(y)/len(y)
    ssreg = sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = sum([ (yi - ybar)**2 for yi in y])
    rsquared = ssreg / sstot
    return (coeffs[0],coeffs[1],rsquared)


if __name__ == """__main__""":
    from numpy import *
    from matplotlib import plot, title, show , legend


    y0=0
    n = 100
    m = 0.75
    err = 10.0
    
    x = linspace(0,100,n)
    y = m*x+random.normal(0,err,n)


    
    print "fixed: ", 
    print linear_reg_fixed_intercept(x,y,y0)

    print "regular " ,
    print linear_reg(x,y)
    
#matplotlib ploting
    title('Linear Regression Example')
    plot(x,y,'g.--')
#    plot(y,xn,'k.')
    plot(x,y,'r.-')
    legend(['original', 'regression'])
    
    show()
