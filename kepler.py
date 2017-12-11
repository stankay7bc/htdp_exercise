#!/usr/bin/python3.5
from math import isclose, sqrt

def area_trapezoid(right,left,afun):
    """ Number Number (Number->Number) -> Number
        compute an area of a trapezoid
    """
    return 0.5*(right-left)*(afun(left)+afun(right))

def integrate_kepler(afun,pa,pb):
    """ (Number->Number) Number Number -> Number
        compute approximation of integral of afun from pa to pb 
      
      EXAMPLES: 
      >>> EPS=0.1
      >>> isclose(integrate_kepler(lambda x: 20,12,22),200,rel_tol=EPS)
      True
      >>> isclose(integrate_kepler(lambda x: 2*x,0,10),100,rel_tol=EPS)
      True
      >>> isclose(integrate_kepler(lambda x: 3*sqrt(x),0,10),1000,rel_tol=EPS)
      False
      >>> isclose(integrate_kepler( 
      ...     lambda x: 3*sqrt(x),0,10),2*(10**(3/2)-0),rel_tol=EPS)
      True
    """
    mid = (pa+pb)/2
    return area_trapezoid(mid,pa,afun)+area_trapezoid(pb,mid,afun)
  

if __name__=="__main__":
    import doctest
    doctest.testmod()
