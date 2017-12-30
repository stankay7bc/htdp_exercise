from functools import reduce
# https://docs.python.org/3.5/library/typing.html
from typing import List, Tuple
from math import isclose

Equation = Tuple[List[float],float]
# represents an linear equation as a tuple 
# that contains a list of coefficients 
# and a right hand side value

SystemOfEquations = List[Equation]
# represents a system of linear equations
# Example: [([1, 1], 2), ([2, 3], 1)]

Solution = List[float]
# represents a solution to a system of linear equations

TriangularSOE = SystemOfEquations
# is a SystemOfEquations
# such that the Equations left sides are decreasing in length:
# n,n-1,...,1
# This data represents triangular matrix

def check_solution(asoe: SystemOfEquations, sol: Solution) -> bool:
    """ 
    return true if sol is a solution for asoe;
    ASSUME that asoe and sol are of the same length.
    
    NOTE: add more tests
    
    EXAMPLES:
    >>> check_solution(
    ...     [([2, 2, 3], 10), 
    ...      ([2, 5, 12], 31), 
    ...      ([4, 1, -2], 1)], [1, 1, 2]) # doctest: +SKIP
    True
    >>> check_solution(
    ...     [([1, 1], 6), 
    ...      ([-3, 1], 2)], [1, 5]) # doctest: +SKIP
    True
    """
    def check_equation(equation: Equation) -> bool:
        # check the equality of the left-hand side of the equation
        # and the right-hand side
        return isclose(
            sum(map(lambda co,va: co*va,equation[0],sol)),
            equation[1],
            abs_tol=1e-6)
        
    return reduce(
        lambda base,equation: base and check_equation(equation),asoe,True) 


def subtract(eq1: Equation, eq2: Equation) -> Equation:
    """ subtract the multiple of eq1 from eq2 to get
        0 in the first position of eq2 and return  
        resulted equation; leave off 0 at first position 
        
    EXAMPLES:
    >>> subtract(([3, 9], 21), ([-3, -8], -19)) # doctest: +SKIP
    ([1.0], 2.0)
    >>> subtract(([2, 3, 3], 8), ([2, 3, -2], 3)) # doctest: +SKIP
    ([0.0, -5.0], -5.0)
    >>> subtract(([0, -5], -5), ([-8, -4], -12)) # doctest: +SKIP
    ([-4.0], -12.0)
    """
    mult = -eq2[0][0]/eq1[0][0] 
    return (
        list(map(lambda el1, el2: mult*el1+el2, eq1[0], eq2[0]))[1:],
        mult*eq1[1]+eq2[1])

def triangulate(soe: SystemOfEquations) -> TriangularSOE:
    """ transform soe into triangular matrix
    
        EXAMPLES:
        >>> soe1 = [([2, 2, 3], 10), 
        ...         ([2, 5, 12], 31), 
        ...         ([4, 1, -2], 1)] 
        >>> soe2 = [([2, 3, 3], 8), 
        ...         ([2, 3, -2], 3), 
        ...         ([4, -2, 2], 4)] 
        >>> triangulate(soe1) # doctest: +NORMALIZE_WHITESPACE +SKIP
        [([2, 2, 3], 10), ([3.0, 9.0], 21.0), ([1.0], 2.0)] 
        >>> triangulate(soe2) # doctest: +NORMALIZE_WHITESPACE +SKIP
        [([2, 3, 3], 8), ([-8.0, -4.0], -12.0), ([-5.0], -5.0)] 
        
    """
    if(len(soe)==0):
        return []
    else:
        allZeros = reduce(lambda base, elem: base and elem[0][0]==0, soe, True)
        if allZeros:
            # report an error
            pass
        else:
            soe = rotate_to_nonzero(soe)
            result = [soe[0]]
            result.extend(triangulate([subtract(soe[0],eq) for eq in soe[1:]]))
            return result

def rotate_to_nonzero(soe: SystemOfEquations) -> SystemOfEquations:
    """ place an equation with first non-zero coefficient at the front;
        ASSUMES the soe contains at least one such equation
    """
    if soe[0][0][0] == 0: 
        return rotate_to_nonzero(soe[1:]+soe[:1])
    else:
        return soe

def solve(soe: TriangularSOE) -> Solution:
    """ find solution for a triangular soe
    
        EXAMPLES:
        >>> tsoe1 = [([2, 3, 3], 8), ([-8.0, -4.0], -12.0), ([-5.0], -5.0)]
        >>> tsoe2 = [([1, 2, 3], 3), ([1.0, 2.0], 9.0), ([-1.0], -11.0)]
        >>> solve(tsoe1) == [1, 1, 1]
        True
        >>> solve(tsoe2) == [-4, -13, 11]
        True
    """

    def walk_soe(soe: TriangularSOE, acc: List) -> Solution:
        if len(soe)==0:
            return acc
        else:
            lhs = sum(map(lambda e1, e2: e1*e2, soe[0][0][1:], acc))
            acc.insert(0,(soe[0][1]-lhs)/soe[0][0][0])
            return walk_soe(soe[1:],acc)

    return walk_soe(list(reversed(soe)),[])
    
def gauss(soe: SystemOfEquations) -> Solution:
    """ solve soe 
    
        EXAMPLES: 
        >>> soe3 = [
        ... ([1, 1, -3], 2),
        ... ([3, -2, 1], -1),
        ... ([2, 1, -2], 0)]
        >>> check_solution(soe3,gauss(soe3))
        True
        
        REFS: https://goo.gl/7NxJfv
    """
    return solve(triangulate(soe))

if __name__=="__main__":
    import doctest
    doctest.testmod()