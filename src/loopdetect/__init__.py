"""
LoopDetect

A Python library to detect feedback loops in ordinary differential equation models.
"""

__version__ = "0.2.0"
__author__ = 'Katharina Baum'
__contributors__ = ["DILiS Lab", "Hang Mai Anh Vo (AG Wolf, MDC)"]

from .core import (
    find_loops,
    find_loops_noscc,
    find_loops_vset,
    find_edge,
    loop_summary,
    compare_loops,
    sort_loop_index,
)

__all__ = [
    "find_loops",
    "find_loops_noscc",
    "find_loops_vset",
    "find_edge",
    "loop_summary",
    "compare_loops",
    "sort_loop_index",
]