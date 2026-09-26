# Assignment 3: Algorithm Efficiency and Scalability

Implementation, analysis, and empirical evaluation of **Randomized Quicksort** and a **Hash Table with Chaining**.
The full write-up (proofs, tables, charts, and discussion) is in **[REPORT.md](REPORT.md)**.

## Repository contents

| File | Purpose |
|---|---|
| `quicksort.py` | Deterministic (first-element pivot), Randomized (2-way), and Randomized 3-way quicksort. All sort in place. |
| `hash_table.py` | `HashTableChaining`: chaining, universal hash family, insert/search/delete, dynamic resizing |
| `test_algorithms.py` | 18 unit tests covering edge cases and correctness against Python's `dict` |
| `benchmark_quicksort.py` | Timing and comparison counts on random, sorted, reverse-sorted, and repeated-value arrays |
| `benchmark_hashing.py` | Load-factor study, scaling with and without resizing, adversarial-key test |
| `REPORT.md` | Theoretical analysis, empirical results, and discussion |
| `results/` | CSV data and PNG charts produced by the benchmarks |

## Requirements

- Python 3.8 or newer
- `matplotlib`, needed only for the charts. The algorithms and tests use the standard library only.

```bash
pip install -r requirements.txt
```

## How to run

```bash
# Quick demos
python quicksort.py
python hash_table.py

# Unit tests
python -m unittest -v

# Benchmarks (write CSV + PNG files to results/)
python benchmark_quicksort.py            # about 2-4 minutes; quadratic cases dominate
python benchmark_quicksort.py --quick    # smoke test, a few seconds
python benchmark_hashing.py              # about 1-2 minutes
python benchmark_hashing.py --quick
```

`benchmark_quicksort.py` supports these options:

- `--repeats N`: number of runs per data point. The median is reported.
- `--budget S`: stop growing n for an algorithm/input pair once a single run takes more than S seconds.
- `--seed`: random seed for input generation.
- `--out DIR`: output directory for results.

### Using the code

```python
from quicksort import randomized_quicksort, randomized_quicksort_3way
data = [5, 3, 8, 3, 1]
randomized_quicksort(data)          # data is now [1, 3, 3, 5, 8]

from hash_table import HashTableChaining
t = HashTableChaining()
t.insert("apple", 1)
t.search("apple")                   # 1
t.delete("apple")                   # True
t["pear"] = 2; "pear" in t          # dict-style access also works
```

## Summary of findings

### Quicksort

- **Randomized quicksort runs in Θ(n log n) expected time on every input with distinct keys.**
  - Measured comparisons matched the exact expectation 2(n+1)Hₙ − 4n. For example, at n = 64,000 the formula predicts 1,234,438 and the measurement was 1,218,612.
  - Sorted, reverse-sorted, and random inputs all gave nearly the same times.
- **Deterministic first-pivot quicksort is Θ(n²) on sorted and reverse-sorted input.**
  - It made exactly n(n−1)/2 comparisons.
  - At n = 16,000 it was about 170× slower than the randomized version, and the gap grows with n.
- **On random input, deterministic quicksort is 10–20% faster.** The randomized version pays for random-number generation, and a random array already gives the first-element pivot a random rank.
- **Pivot randomization does not help with many duplicates when a 2-way partition is used.** Both 2-way versions were quadratic on arrays containing only 10 distinct values. The 3-way variant was about 800× faster on that input and only 20–30% slower on distinct keys.

### Hash table with chaining

- **Measured search costs matched simple-uniform-hashing theory within about 3%.**
  - Unsuccessful searches examined α keys on average.
  - Successful searches examined 1 + α/2 keys on average.
  - Both grow linearly in the load factor α.
- **Dynamic resizing keeps operations near constant time as n grows.** The table doubles at α > 0.75 and halves at α < 0.125.
  - With resizing, α stayed at about 0.5 up to n = 1,000,000, so the expected chain length did not change.
  - A fixed table of 1,024 slots slowed about 18× by n = 256,000.
  - The small rise in time for the resizing table at large n comes from CPU cache effects, not from longer chains.
- **The universal hash function ((a·k + b) mod p) mod m prevents structured keys from causing collisions.** On keys that are all multiples of m:
  - `k mod m` placed all 4,096 keys in one slot.
  - The universal family averaged 7,497 colliding pairs, below its theoretical bound of 8,190.

Timings were measured on a Linux machine with Python 3.12. Absolute numbers will differ on other machines, but growth rates and ratios should be similar.
