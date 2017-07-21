"""Depth functions or smth
"""

from __future__ import division
from positron_constants import p_depth,z
from slice_sampling import slice_sampling
from collections import OrderedDict

class DepthFunctions:
    """A simple lookup table for depth functions"""
    def __init__(self, accuracy, cache_size=10**6):
        self.depths = OrderedDict()
        self.accuracy = float(accuracy)
        self.cache_size = cache_size
        self.element_count = 0
    def __getitem__(self, energy):
        unique_index = int(round(energy / self.accuracy))
        #if energy is 0, it is rounded up to the accuracy value
        if unique_index == 0: unique_index = 1
        #If it's in the cache, stick it at the end and return the next number
        if unique_index in self.depths:
            gen = self.depths[unique_index]
            del self.depths[unique_index]
            self.depths[unique_index] = gen
            return next(gen)
        #If our cache isn't full, add it
        elif self.element_count < self.cache_size:
            gen = slice_sampling(p_depth(unique_index*self.accuracy),z,width=10000)
            self.depths[unique_index] = gen
            self.element_count += 1
            return next(gen)
        #If our cache is full and our wanted generator is not in it,
        # remove the generator in the cache that was least recently used 
        # and add the new value
        else:
            gen = slice_sampling(p_depth(unique_index*self.accuracy),z,width=10000)
            self.depths.popitem()
            self.depths[unique_index] = gen
            return next(gen)