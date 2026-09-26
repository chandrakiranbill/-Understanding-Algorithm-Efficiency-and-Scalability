"""
Quicksort implementations for Assignment 3.

Three variants are provided so that the effect of *pivot choice* and the effect
of *partition scheme* can be studied separately:

    deterministic_quicksort   -- pivot = first element, 2-way (Lomuto) partition
    randomized_quicksort      -- pivot chosen uniformly at random, 2-way partition
    randomized_quicksort_3way -- pivot chosen uniformly at random, 3-way partition
                                 (recommended; linear-time on arrays of equal keys)

All functions sort the list IN PLACE and also return it for convenience.
Every variant recurses only into the smaller side of a partition and loops on
the larger side, so the recursion depth is O(log n) even when the running time
degrades to O(n^2). This means no variant can crash with a RecursionError on
already-sorted input, which lets us measure worst-case behaviour safely.

Each function accepts an optional `counter` dict. If given, the number of
element comparisons against the pivot is accumulated in counter["comparisons"].
Comparisons are a machine-independent cost measure that can be compared
directly with the theoretical bound of about 2n ln n.
"""

import random

__all__ = [
    "deterministic_quicksort",
    "randomized_quicksort",
    "randomized_quicksort_3way",
]


# ---------------------------------------------------------------------------
# Partition routines
# ---------------------------------------------------------------------------
def _lomuto_partition(a, lo, hi, counter):
    """Partition a[lo..hi] around the pivot stored at a[lo].

    Returns the final index of the pivot. Elements <= pivot end up on its
    left, elements > pivot on its right.
    """
    pivot = a[lo]
    i = lo  # a[lo+1..i] holds the elements <= pivot
    for j in range(lo + 1, hi + 1):
        if a[j] <= pivot:
            i += 1
            a[i], a[j] = a[j], a[i]
    a[lo], a[i] = a[i], a[lo]
    if counter is not None:
        counter["comparisons"] = counter.get("comparisons", 0) + (hi - lo)
    return i


def _three_way_partition(a, lo, hi, counter):
    """Dijkstra's Dutch-national-flag partition of a[lo..hi] around a[lo].

    Returns (lt, gt) such that
        a[lo..lt-1] < pivot,  a[lt..gt] == pivot,  a[gt+1..hi] > pivot.
    """
    pivot = a[lo]
    lt, i, gt = lo, lo + 1, hi
    comps = 0
    while i <= gt:
        x = a[i]
        comps += 1
        if x < pivot:
            a[lt], a[i] = x, a[lt]
            lt += 1
            i += 1
        elif x > pivot:
            a[i], a[gt] = a[gt], x
            gt -= 1
        else:
            i += 1
    if counter is not None:
        counter["comparisons"] = counter.get("comparisons", 0) + comps
    return lt, gt


# ---------------------------------------------------------------------------
# Public sorting functions
# ---------------------------------------------------------------------------
def deterministic_quicksort(a, counter=None):
    """Quicksort with the FIRST element of each subarray as the pivot."""
    if a is None:
        raise TypeError("input must be a list, not None")
    lo, hi = 0, len(a) - 1
    stack = []
    while True:
        while lo < hi:
            p = _lomuto_partition(a, lo, hi, counter)
            # Push the larger side, continue with the smaller side.
            if p - lo < hi - p:
                stack.append((p + 1, hi))
                hi = p - 1
            else:
                stack.append((lo, p - 1))
                lo = p + 1
        if not stack:
            return a
        lo, hi = stack.pop()


def randomized_quicksort(a, counter=None, rng=random):
    """Quicksort with a pivot chosen uniformly at random from the subarray.

    Uses the same 2-way partition as `deterministic_quicksort`; the ONLY
    difference is how the pivot is chosen.
    """
    if a is None:
        raise TypeError("input must be a list, not None")
    lo, hi = 0, len(a) - 1
    stack = []
    randint = rng.randint
    while True:
        while lo < hi:
            r = randint(lo, hi)          # uniform random pivot index
            a[lo], a[r] = a[r], a[lo]    # move pivot to the front
            p = _lomuto_partition(a, lo, hi, counter)
            if p - lo < hi - p:
                stack.append((p + 1, hi))
                hi = p - 1
            else:
                stack.append((lo, p - 1))
                lo = p + 1
        if not stack:
            return a
        lo, hi = stack.pop()


def randomized_quicksort_3way(a, counter=None, rng=random):
    """Randomized quicksort with 3-way partitioning.

    Keys equal to the pivot are placed in the middle and never looked at
    again, so arrays with many repeated values sort in O(n log k) expected
    time, where k is the number of distinct keys (O(n) if all keys are equal).
    """
    if a is None:
        raise TypeError("input must be a list, not None")
    lo, hi = 0, len(a) - 1
    stack = []
    randint = rng.randint
    while True:
        while lo < hi:
            r = randint(lo, hi)
            a[lo], a[r] = a[r], a[lo]
            lt, gt = _three_way_partition(a, lo, hi, counter)
            left_size, right_size = lt - lo, hi - gt
            if left_size < right_size:
                stack.append((gt + 1, hi))
                hi = lt - 1
            else:
                stack.append((lo, lt - 1))
                lo = gt + 1
        if not stack:
            return a
        lo, hi = stack.pop()


if __name__ == "__main__":
    demo = [5, 3, 8, 3, 1, 9, 2, 3, 7, 0]
    print("input:                ", demo)
    print("deterministic:        ", deterministic_quicksort(demo[:]))
    print("randomized:           ", randomized_quicksort(demo[:]))
    print("randomized (3-way):   ", randomized_quicksort_3way(demo[:]))
