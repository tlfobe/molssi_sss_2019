"""
geometry analysis functions
"""

import os
import numpy as np

def calculate_distance(rA, rB):
    """Calculate the distance between points A and B. Assumes rA and rB are numpy arrays."""
    dist_vec = (rA-rB)
    distance = np.linalg.norm(dist_vec)
    return distance

def calculate_distance_list(rA, rB):
    """
    Calculate the distance between points A and B. Assums rA and rB are lists.
    
    Iterates over the dimensions of both rA and rB summing over their squares of each element of A - B and takes the square root of their sum.
    
    Parameters
    ----------
    rA : array_like
        Array of xyz coordinates of point A
    rB : array_like
        Array of xyz coordinate of point B

    Returns
    -------
    distance : float
        Distance between points A and B

    Examples
    --------
    rA = np.array([1,0,0])
    rB = np.array([0,0,0])
    dist = calculate_distance_list(rA, rB)

    print(dist)
    >>> 1.0
    """
    squared_sum = 0
    for dim in range(len(rA)):
        squared_sum += (rA[dim] - rB[dim])**2
    
    distance = np.sqrt(squared_sum)
    return distance


def build_bond_list(coordinates, max_bond=1.5, min_bond=0):
    """
    Build list of bonsd from atomic coordinates based on distance.

    Parameters
    ----------

    coodinates : np.array
        An Array of atomic coordinates. Size should be (n,3) where n is the number of particles
    max_bond : float, optional
        The maximum distance between atoms to be considered a bond. Default is 2.93 bohr
    min_bond : float, optional
        The minimum distance between atoms to be considered a bond.

    Returns
    -------

    bonds : dict
        A dictionary of tuples of atom indices which points to a bond distance

    Examples
    --------

    >>> coords = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]])
    >>> build_bond_list(coords)
    {(0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0, (1, 2): 1.4142135623730951, (1, 3): 1.4142135623730951, (2, 3): 1.4142135623730951}


    """
    num_atoms = len(coordinates)
    
    bonds = {}
    
    for atom1 in range(num_atoms):
        for atom2 in range(atom1, num_atoms):
            distance = calculate_distance(coordinates[atom1], coordinates[atom2])
            
            if distance > min_bond and distance < max_bond:
                bonds[(atom1, atom2)] = distance 
    
    return bonds


