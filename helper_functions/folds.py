"""A module for composable left folds over any iterable to allow processing
of files in constant space with a single pass.
This is very unpythonic, but it'll work. Inspired by the folds library
for Haskell by Gabriel 'Tekmo' Gonzalez

Copyright 2016 Robert 'Probie' Offner
This is free for anyone to use and modify without attribution if you've got
nothing better to do or something. Maybe whilst you're at it you could fix
python2 so it could memmap files larger than 2gb or something.
"""

from __future__ import division
from functools import partial
from itertools import imap
from collections import Counter
import operator
import math
import numpy as np

class Fold:
    """A simple left fold. Do not use directly, but use a FoldBuilder instead
       since a Fold is only good for one use
    """

    def __init__(self, f, seed, done):
        self.f = f
        self.seed = seed
        self.done = done

    def __call__(self,xs):
        return self.done(reduce(self.f,xs,self.seed))

class FoldBuilder:
    """A class for building Folds.
    """
    def __init__ (self, f, seed, done=None):
        self.f = f
        self.seed = seed
        self.done = _identity if done is None else done

    def __call__(self,xs):
        """This builds a fold, and then uses the fold's __call__ method on
        the argument given
        """
        
        return Fold(self.f,self.seed,self.done)(xs)

    def __and__(self,other):
        """Compose two folds together into a single fold
        """
        def combine_done((x,y)):
            done_function = self.done(x)
            try:
                return done_function(other.done(y))
            except TypeError:
                # We probably needed partial application
                # There's no general way to check this, and I'm sorry
                # but parsing docstrings is not the solution. Plus just trying
                # and letting it fail is pythonic, right?
                return partial(done_function,other.done(y))
        return FoldBuilder(lambda (s1,s2), x: (self.f(s1,x), other.f(s2,x)),
                           (self.seed,other.seed),
                           done=combine_done)

    def __rand__(self,other):
        """Compose a function and a fold together into a single fold
        """
        if not callable(other):
            raise TypeError ('unsupported operand type(s) for &: \'' +
                             type(other).__name__ + '\' and \'fold\'')
        return lift_function_into_fold(other) & self
    def __rshift__(self,other):
        """Combine two folds, but only return the results of the second.
        This is mainly useful when the first fold is side-effecting
        """
        return FoldBuilder(lambda (s1,s2), x:(self.f(s1,x), other.f(s2,x)),
                           (self.seed, other.seed),
                           done=lambda (_,y):other.done(y))

def _identity(x):
    return x

def _binary_no_op(x,y):
    return None
    
def _const(x):
    """Return a constant value"""
    return lambda _: x

def lift_function_into_fold(f):
    return FoldBuilder(_binary_no_op, None, done=_const(f))

# Let's just define a few simple folds for people to use

length = FoldBuilder(lambda x, _: 1+x,0)
sum_generator = FoldBuilder(operator.add,0)
average = (operator.truediv) & sum_generator & length

minimum = FoldBuilder(lambda x,y: y if x is None else min(x,y), None)
maximum = FoldBuilder(max,None)

def maketuple(x,y):
    return (x,y)

maxandmin = maketuple & maximum & minimum

def when(p,fold):
    """ Fold only over the elements for which p
    holds
    """
    f = fold.f
    z = fold.seed
    done = fold.done
    def worker(acc,x):
        if p(x): return f(acc,x)
        return acc
    return FoldBuilder(worker,z,done=done)

def before(fold,g):
    """ Apply the function 'g' before each element
    of the fold
    """

    f = fold.f
    z = fold.seed
    done = fold.done
    def worker(acc,x):
        return f(acc,g(x))
    return FoldBuilder(worker,z,done=done)


def print_function(val):
    print val

def side_effect_fold(f, done=None):
    return FoldBuilder(lambda _,val: f(val), None, done=done)

printfold = side_effect_fold(print_function)

def collect_into(container):
    return side_effect_fold(container.append)

def collect_list():
    x = []
    return side_effect_fold(x.append,done=_const(x))

def collect_numpy(size):
    out = np.empty(size)
    iterator = np.nditer(out, op_flags=['writeonly'])
    def worker(val):
        try:
            next(iterator)[...] = val
        except StopIteration:
            raise ValueError('Generator larger than numpy array in fold')
    return side_effect_fold(worker, done=lambda _: out)

def _tuple_combine(acc, f):
    return maketuple & acc & f

def _addWithDepth(n,xs,result):
    if n == 1:
        result.append(xs)
        return
    _addWithDepth(n-1,xs[0],result)
    result.append(xs[1])

def select(*args):
    def selection_function(item):
        if len(args) == 1: return item[args[0]]
        return tuple(map(partial(operator.getitem,item),args))
    return selection_function

def combine(*args):
    """Combine the given folds into a single fold returning a list"""
    new_fold = reduce(_tuple_combine,args)
    def combine_done(old_done):
        def worker(fold_result):
            result = []
            result_tuple = old_done(fold_result)
            _addWithDepth(len(args),result_tuple,result)
            return result
        return worker
    new_fold.done = combine_done(new_fold.done)
    return new_fold

def tally():
    """Return a fold builder which tallies the elements"""
    res = Counter()
    def inc(x):
        res[x] +=1
    return side_effect_fold(inc, done=_const(res))

def histogram(bin_count, min_value, max_value, normalised=False, bin_dist=None):
    """Make a histogram with bin_count evenly spaced bins between
    min_value and max_value. Anything less than min_value will
    go in the smallest bin, and anything greater than max_value
    will be in the largest bin. By default the values in the bins
    is the number of occurences. To change this behaviour, set
    normalised to True
    """
    def map_to_bin_number(val):
        val = max(min_value,min(max_value,val)) # clamp
        return int(round((val-min_value) / (max_value - min_value)
                         * (bin_count - 1)))
    
    if normalised:
        return partial(maketuple,np.linspace(min_value,max_value,bin_count)) \
            & ((lambda x, y: np.array(map(lambda a: a / x,y))) & length &
               ((lambda x: ([x[n] for n in xrange(bin_count)])) &
                before(tally(),map_to_bin_number)))
    return partial(maketuple,np.linspace(min_value,max_value,bin_count)) \
        & ((lambda x: np.array([x[n] for n in xrange(bin_count)]))&
           before(tally(),map_to_bin_number))
