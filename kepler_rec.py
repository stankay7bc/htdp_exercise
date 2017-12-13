#!/usr/bin/python3.5
from math import isclose, sqrt
from kepler import integrate_kepler, area_trapezoid

BASIC_INT = 0.005 

def integrate_dc(afun,pa,pb):
    """ (Number -> Number) Number Number -> Number
    compute integral of afun from a to b
    by recursively applying kepler technique
    
    EXAMPLES:
    >>> EPS=0.001
    >>> isclose(integrate_dc(lambda x: 20,12,22),200,rel_tol=EPS)
    True
    >>> isclose(integrate_dc(lambda x: 2*x,0,10),100,rel_tol=EPS)
    True
    >>> isclose(integrate_dc(lambda x: 3*sqrt(x),0,10),1000,rel_tol=EPS)
    False
    >>> isclose(integrate_dc( 
    ...     lambda x: 3*sqrt(x),0,10),2*(10**(3/2)-0),rel_tol=EPS)
    True
    """
    if (pb-pa) <= BASIC_INT:
        return integrate_kepler(afun,pa,pb)
    else:
        mid = (pb+pa)/2
        return integrate_dc(afun,pa,mid) + \
            integrate_dc(afun,mid,pb)

def integrate_adaptive(afun,pa,pb,eps=0.001):
    """ adaptive version of integrate_dc
    
    EXAMPLES:
    >>> EPS=0.001
    >>> isclose(integrate_adaptive(lambda x: 20,12,22),200,rel_tol=EPS)
    True
    >>> isclose(integrate_adaptive(lambda x: 2*x,0,10),100,rel_tol=EPS)
    True
    >>> isclose(integrate_adaptive(lambda x: 3*sqrt(x),0,10),1000,rel_tol=EPS)
    False
    >>> isclose(integrate_adaptive( 
    ...     lambda x: 3*sqrt(x),0,10),2*(10**(3/2)-0),rel_tol=EPS)
    True
    """
    appr1 = area_trapezoid(pb,pa,afun)
    appr2 = integrate_kepler(afun,pa,pb)
    if abs(appr1-appr2)<=eps*(pb-pa):
        return appr2
    else:
        mid = (pb+pa)/2
        return integrate_adaptive(afun,pa,mid) + \
            integrate_adaptive(afun,mid,pb)


if __name__=="__main__":
    import doctest
    doctest.testmod()