#!/usr/bin/python3.5
from math import isclose, sqrt

def integrateRS(afun,a,b,R=10e0):
    """ (Number -> Number) Number Number -> Number
        use Riemann sum to compute 
        integral of afun from a to b
    EXAMPLES:
    >>> EPS=0.01
    >>> isclose(integrateRS(lambda x: 20,12,22),200,rel_tol=EPS)
    True
    >>> isclose(integrateRS(lambda x: 2*x,0,10),100,rel_tol=EPS)
    True
    >>> isclose(integrateRS(lambda x: 3*sqrt(x),0,10),1000,rel_tol=EPS)
    False
    >>> isclose(integrateRS( 
    ...     lambda x: 3*sqrt(x),0,10),2*(10**(3/2)-0),rel_tol=EPS)
    True
    """
    width = (b-a)/R
    return width*sum(afun(a+width/2+width*i) for i in range(0,int(R)))

if __name__=="__main__":
    import doctest
    doctest.testmod()