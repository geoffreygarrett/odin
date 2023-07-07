# import pytest
import numpy as np
import pydin
import os
from llvmlite import binding
import numba

# load the library into LLVM
path = os.path.abspath(pydin.core.__file__)
binding.load_library_permanently(path)
from numba.experimental import jitclass
import numpy as np

state_test = (
    np.array([0, 1], dtype=np.int32),
    np.array([0.5, 0.75, 200.0], dtype=np.float64)
)

from numba import types, int32, float64
# from numba.types import Tuple
from typing import Tuple

# Specifying the data layout for the jitclass
spec = [
    ('raw_state', types.Tuple([types.Array(int32, 1, 'C'), types.Array(float64, 1, 'C')])),
]


@jitclass(spec)
class MIPState:
    def __init__(self, raw_state: types.Tuple([types.Array(int32, 1, 'C'), types.Array(float64, 1, 'C')])):
        self.raw_state = raw_state

    @staticmethod
    def from_state(raw_state: types.Tuple([types.Array(int32, 1, 'C'), types.Array(float64, 1, 'C')])):
        return MIPState(raw_state)

    def integer(self):
        return self.raw_state[0]

    def continuous(self):
        return self.raw_state[1]

    def to_raw(self):
        return self.raw_state


spec2 = [
    ('mip_state', MIPState.class_type.instance_type),
]

from numba.types import Tuple


@jitclass(spec2)
class Step1State:
    def __init__(self, raw_state: types.Tuple([types.Array(int32, 1, 'C'), types.Array(float64, 1, 'C')])):
        self.mip_state = MIPState(raw_state)

    def integer(self):
        return self.mip_state.integer()

    def continuous(self):
        return self.mip_state.continuous()

    def to_raw(self):
        return self.mip_state.to_raw()

    @property
    def epoch_mjd(self):
        return self.mip_state.continuous()[0]


state = Step1State(
    (
        np.array([0, 1], dtype=np.int32),
        np.array([0.5, 0.75, 200.0], dtype=np.float64)
    ))


@numba.njit()
def test_python1(state: types.Tuple([types.Array(int32, 1, 'C'), types.Array(float64, 1, 'C')])):
    for i in range(state[0].size):
        print("integer[", i, "]:", state[0][i])
    for i in range(state[1].size):
        print("continuous[", i, "]:", state[1][i])
    print("int_size: ", state[0].size)
    print("continuous_size: ", state[1].size)


#
# @numba.njit()
# def test_python2(state: Tuple[np.ndarray[np.int32], np.ndarray[np.float64]]):
#     m_state = Step1State(state)
#     for i in range(m_state.integer().size):
#         print("integer[", i, "]:", m_state.integer()[i])
#     for i in range(m_state.continuous().size):
#         print("continuous[", i, "]:", m_state.continuous()[i])
#     print("int_size: ", m_state.integer().size)
#     print("continuous_size: ", m_state.continuous().size)
#
#
import time


#

def time_test_performance(func):
    start_time = time.time()
    func()
    end_time = time.time()
    return end_time - start_time


#
# @numba.njit()
# def test_performance1():
#     pass
#
#
# @numba.njit()
# def test_performance2():
#     for i in range(10000000):
#         state_test = (
#             (0, 1),
#             (0.5, 0.75, 200.0)
#         )
#         m_state = Step1State(state_test)
#         x = m_state.epoch_mjd - 1.0
#
#
# @numba.njit()
# def test_performance3():
#     for i in range(10000000):
#         state_test = (
#             (0, 1),
#             (0.5, 0.75, 200.0)
#         )
#         x = state_test[0][0] - 1.0
#
#
# print("test_performance")
# print(time_test_performance(test_performance1))
# print(time_test_performance(test_performance2))
# print(time_test_performance(test_performance3))
# print()

test_python1(state_test)
pydin.test_cpp(state_test)
pydin.test_cpp(state.to_raw())

print(state.epoch_mjd)
import timeit
