'''
 write down a Jacobian method.
   Then compute the Jacobian matrix of these functions wrt the variables. 
   Write a rank computing method. 
   The inputs for these two methods should be a list of functions and a list of variables.
     Next write an evaluation method, which takes a list of functions, and a list of floating numbers. 
     This method computes the Jacobian with these inputs, and compute the rank.
'''

import sympy as sp
from sympy import Matrix, symbols, diff, lambdify, simplify
import numpy as np
from typing import List, Union, Tuple
import warnings

def compute_jacobian(functions: List[sp.Expr], 
                     variables: List[sp.Symbol]) -> sp.Matrix:
    """
    Compute the Jacobian matrix of a list of functions with respect to variables.
    
    The Jacobian matrix J is defined as:
    J[i,j] = df_i/dx_j
    
    Parameters
    ----------
    functions : List[sp.Expr]
        List of functions f_1, f_2, ..., f_m
    variables : List[sp.Symbol]
        List of variables x_1, x_2, ..., x_n
        
    Returns
    -------
    sp.Matrix
        Jacobian matrix of shape (m × n) where m = len(functions), n = len(variables)
        
    """
    m = len(functions)
    n = len(variables)
    
    # Initialize Jacobian matrix
    J = sp.zeros(m, n)
    
    # Compute Jacobian
    for i, func in enumerate(functions):
        for j, var in enumerate(variables):
            J[i, j] = diff(func, var)
    
    return J


def compute_rank(jacobian: sp.Matrix, 
                         simplify_result: bool = True) -> int:
    """
    Compute the rank Jacobian.
    
    Parameters
    ----------
    jacobian : sp.Matrix
        Jacobian matrix
        
    Returns
    -------
    int
        Rank of the Jacobian matrix
    """
    if simplify_result:
        J_simplified = jacobian.applyfunc(simplify)
    else:
        J_simplified = jacobian
    
    return J_simplified.rank()



def evaluate_jacobian_at_point(functions: List[sp.Expr],
                               variables: List[sp.Symbol],
                               values: List[Union[float, complex]],
                               tol: float = 1e-10) -> Tuple[np.ndarray, int, List[float]]:
    """
    Evaluate Jacobian at a specific point and compute its rank.
    
    Parameters
    ----------
    functions : List[sp.Expr]
        List of functions
    variables : List[sp.Symbol]
        List of variables
    values : List[Union[float, complex]]
        List of numerical values to substitute for variables
    tol : float, optional
        Tolerance for rank computation (default: 1e-10)
        
    Returns
    -------
    jacobian_numeric : np.ndarray
        Evaluated Jacobian matrix as numpy array
    rank : int
        Numerical rank of the Jacobian
    singular_values : List[float]
        List of singular values (for diagnostic purposes)
        
    Raises
    ------
    ValueError
        If length of values doesn't match length of variables
    """
    if len(values) != len(variables):
        raise ValueError(f"Number of values ({len(values)}) must match "
                        f"number of variables ({len(variables)})")
    
    # Compute Jacobian
    print("Computing symbolic Jacobian...")
    J_symbolic = compute_jacobian(functions, variables)
    
    # Create substitution
    subs_dict = dict(zip(variables, values))
    
    # Evaluate Jacobian at the point
    print("Evaluating Jacobian at specified point...")
    J_evaluated = J_symbolic.subs(subs_dict)
    
    # Convert to numpy array (handling complex numbers)
    print("Converting to numeric array...")
    J_numeric = np.array(J_evaluated.tolist(), dtype=complex)
    
    # Compute rank 
    print("Computing rank")
    singular_values = np.linalg.svd(J_numeric, compute_uv=False)
    rank = np.sum(singular_values > tol)
    
    return J_numeric, rank, singular_values.tolist()


def batch_evaluate_jacobian(functions: List[sp.Expr],
                           variables: List[sp.Symbol],
                           value_sets: List[List[Union[float, complex]]],
                           tol: float = 1e-10) -> List[dict]:
    """
    Evaluate Jacobian rank at multiple points.
    
    Useful for statistical analysis of rank across parameter space.
    
    Parameters
    ----------
    functions : List[sp.Expr]
        List of symbolic functions
    variables : List[sp.Symbol]
        List of symbolic variables
    value_sets : List[List[Union[float, complex]]]
        List of parameter value sets to evaluate
    tol : float, optional
        Tolerance for rank computation (default: 1e-10)
        
    Returns
    -------
    List[dict]
        List of dictionaries containing results for each evaluation:
        - 'values': parameter values used
        - 'rank': computed rank
        - 'singular_values': list of singular values
        - 'condition_number': condition number of Jacobian
        
    Examples
    --------
    >>> x, y = symbols('x y')
    >>> f1 = x**2 + y**2
    >>> f2 = x*y
    >>> values = [[1.0, 2.0], [0.5, 1.5], [2.0, 3.0]]
    >>> results = batch_evaluate_jacobian([f1, f2], [x, y], values)
    """
    # Compute symbolic Jacobian once
    print("Computing symbolic Jacobian...")
    J_symbolic = compute_jacobian(functions, variables)
    
    # Convert to numerical function using lambdify for speed
    print("Creating numerical function...")
    J_func = lambdify(variables, J_symbolic, modules=['numpy'])
    
    results = []
    
    for i, values in enumerate(value_sets):
        print(f"\nEvaluating point {i+1}/{len(value_sets)}...")
        
        try:
            # Evaluate Jacobian
            J_numeric = np.array(J_func(*values), dtype=complex)
            
            # Compute singular values
            singular_values = np.linalg.svd(J_numeric, compute_uv=False)
            
            # Compute rank
            rank = np.sum(singular_values > tol)
            
            # Compute condition number
            if singular_values[-1] > tol:
                condition_number = singular_values[0] / singular_values[-1]
            else:
                condition_number = np.inf
            
            results.append({
                'values': values,
                'rank': int(rank),
                'singular_values': singular_values.tolist(),
                'condition_number': float(condition_number)
            })
            
        except Exception as e:
            warnings.warn(f"Error evaluating point {i+1}: {e}")
            results.append({
                'values': values,
                'rank': None,
                'singular_values': None,
                'condition_number': None,
                'error': str(e)
            })
    
    return results
