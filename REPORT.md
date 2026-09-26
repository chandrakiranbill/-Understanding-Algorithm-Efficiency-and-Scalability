# Assignment 3 Report: Understanding Algorithm Efficiency and Scalability

**Algorithms studied:** Randomized Quicksort and Hashing with Chaining
**Language:** Python 3.12 (standard library only; `matplotlib` for charts)

All numbers in this report come from `benchmark_quicksort.py` and `benchmark_hashing.py`, and the raw data is in `results/`. Absolute times depend on the machine, so rerunning the scripts will give different milliseconds. The growth rates and ratios discussed below should be the same on any machine.

---

## Part 1: Randomized Quicksort

### 1.1 Implementation

`quicksort.py` provides three in-place variants. They share as much code as possible, so each comparison isolates one design decision.

| Function | Pivot rule | Partition scheme |
|---|---|---|
| `deterministic_quicksort` | first element of the subarray | 2-way (Lomuto) |
| `randomized_quicksort` | uniformly random element of the subarray | 2-way (Lomuto) |
| `randomized_quicksort_3way` | uniformly random element of the subarray | 3-way (Dijkstra "Dutch national flag") |

The deterministic and randomized 2-way versions differ **only** in pivot selection. The randomized version draws an index `r` uniformly from `[lo, hi]`, swaps `a[r]` to the front, and then calls the same partition routine. Any difference between them is therefore caused by the pivot rule alone.

**Design decisions and edge cases**

- **Bounded recursion depth.** Every variant processes the smaller side of each partition first and handles the larger side with an explicit loop and stack. Because the smaller side has at most half the elements, the pending work is O(log n) deep even when running time degrades to Θ(n²). This matters in Python, where the default recursion limit is about 1,000. A naive recursive deterministic quicksort crashes with `RecursionError` on a sorted array of only a few thousand elements. Here the quadratic cases can be measured safely.
- **Empty and single-element arrays.** The loop condition `lo < hi` is false immediately, so the list is returned unchanged.
- **Already sorted and reverse-sorted arrays.** These are handled correctly by all variants. The randomized variants avoid the quadratic behaviour (Section 1.3).
- **Repeated elements.** All variants sort them correctly. The 3-way variant also stays fast: keys equal to the pivot are grouped in the middle and never examined again. The 2-way variants degrade here even with a random pivot, as explained in Sections 1.2 and 1.3.
- **Other inputs.** Any mutually comparable values work, including negatives, floats, and strings. Passing `None` raises `TypeError`.
- **Instrumentation.** An optional `counter` dictionary records the number of elements compared against a pivot. This gives a hardware-independent cost that can be checked against the theoretical formula.

`test_algorithms.py` contains 18 unit tests covering all of these cases. All tests pass.

### 1.2 Average-case analysis

Let the input contain n **distinct** keys, and let z₁ < z₂ < … < zₙ be those keys in sorted order. The running time of quicksort is O(n + X), where X is the total number of comparisons made during partitioning. Each partition call does O(1) work plus one comparison per element examined, and there are at most n partition calls. It is therefore enough to compute E[X].

#### Method 1: Indicator random variables

Define the indicator random variable

  X_ij = I{ zᵢ is compared with zⱼ at some point during the algorithm },  for 1 ≤ i < j ≤ n.

Every comparison is between some element and the current pivot. After a partition, the pivot is in its final position and never takes part in another comparison. So any pair is compared at most once, and

  X = Σᵢ₌₁ⁿ⁻¹ Σⱼ₌ᵢ₊₁ⁿ X_ij,  so by linearity of expectation  E[X] = Σᵢ Σⱼ Pr{zᵢ is compared with zⱼ}.

**Computing the probability.** Let Z_ij = {zᵢ, zᵢ₊₁, …, zⱼ}, which contains j − i + 1 keys.

- As long as no pivot has been chosen from Z_ij, all of its keys lie in the same subarray. A pivot outside Z_ij is either smaller than all of them or larger than all of them, so it cannot separate them.
- When the first pivot from Z_ij is chosen, it is chosen uniformly from the subarray that contains all of Z_ij. Every key of Z_ij is therefore equally likely to be that first pivot, with probability 1/(j − i + 1).
- If that first pivot is zᵢ or zⱼ, the two keys are compared, because the pivot is compared with every element of its subarray.
- If it is some z_k with i < k < j, then zᵢ and zⱼ end up on opposite sides of z_k and are never compared.

Hence

  Pr{zᵢ is compared with zⱼ} = 2 / (j − i + 1).

**Summing.** Substitute k = j − i:

  E[X] = Σᵢ₌₁ⁿ⁻¹ Σₖ₌₁ⁿ⁻ⁱ 2/(k + 1) < Σᵢ₌₁ⁿ⁻¹ 2Hₙ = O(n log n),

where Hₙ = 1 + 1/2 + … + 1/n = ln n + O(1) is the n-th harmonic number. For a lower bound, the pairs with i ≤ n/2 each contribute at least 2(H_{n/2} − 1), so E[X] = Ω(n log n) as well. The expected running time is therefore **Θ(n log n)**.

#### Method 2: Recurrence (exact constant)

Let T(n) be the expected number of comparisons on n distinct keys. The pivot's rank q is uniform on {0, …, n−1}, and partitioning uses n − 1 comparisons, so

  T(n) = (n − 1) + (1/n) Σ_{q=0}^{n−1} [T(q) + T(n−1−q)] = (n − 1) + (2/n) Σ_{q=0}^{n−1} T(q),  with T(0) = T(1) = 0.

To solve it, multiply by n, write the same equation for n − 1, and subtract:

  nT(n) − (n−1)T(n−1) = 2(n−1) + 2T(n−1),  which gives  nT(n) = (n+1)T(n−1) + 2(n−1).

Divide both sides by n(n+1), and use 2(n−1)/(n(n+1)) = 4/(n+1) − 2/n:

  T(n)/(n+1) = T(n−1)/n + 4/(n+1) − 2/n.

Telescoping from 1 to n gives T(n)/(n+1) = 2Hₙ + 4/(n+1) − 4, and therefore

  **T(n) = 2(n+1)Hₙ − 4n ≈ 2n ln n ≈ 1.386 · n log₂ n.**

This exact formula is the dashed "theory" line in `results/quicksort_comparisons.png`.

#### Why randomization matters

The expectation above is taken over the algorithm's **own random choices**, not over a distribution of inputs. The bound therefore holds for **every** input of distinct keys, including sorted and reverse-sorted arrays. Deterministic first-pivot quicksort reaches Θ(n log n) only when the input itself is randomly ordered. On a sorted array its pivot is always the minimum, so each partition removes only one element, and the total number of comparisons is

  (n−1) + (n−2) + … + 1 = n(n−1)/2 = Θ(n²).

Randomized quicksort can still be quadratic on some sequence of coin flips, but that happens with vanishingly small probability. There is no longer any "bad input", only extremely unlucky random draws.

#### The distinct-keys assumption

The analysis assumes that all keys are distinct. With a 2-way (Lomuto) partition, every element that is ≤ the pivot goes to the left. When all keys are equal, every partition therefore places all remaining elements on one side, and the cost is n(n−1)/2 **no matter which pivot is chosen**. Randomization cannot help here. A 3-way partition fixes this, because every key equal to the pivot is finished after one pass. With k distinct values its expected cost is O(n log k), and O(n) when all keys are equal.

### 1.3 Empirical comparison

**Setup.** Sizes n = 500 to 64,000, doubling each time. Each cell reports the **median of 3 runs**, and every run uses a fresh input. The four distributions are:

- random integers in [0, 10⁹]
- already sorted: `0, 1, …, n−1`
- reverse sorted
- repeated elements: random values from {0, …, 9}

Once a single run of an (algorithm, distribution) pair took longer than 4 seconds, larger sizes for that pair were skipped. Correctness of every output was asserted. The tables show median time with the number of comparisons in parentheses.

**Random arrays**

| n | Deterministic (first pivot) | Randomized (2-way) | Randomized (3-way) |
|---:|---:|---:|---:|
| 1,000 | 0.8 ms (10,118) | 1.2 ms (10,252) | 1.4 ms (11,294) |
| 4,000 | 5.0 ms (52,386) | 5.6 ms (54,945) | 6.9 ms (55,931) |
| 16,000 | 21.1 ms (258,702) | 25.7 ms (260,615) | 33.1 ms (263,829) |
| 64,000 | 106.9 ms (1,198,354) | 122.2 ms (1,218,612) | 157.4 ms (1,261,820) |

**Already sorted arrays**

| n | Deterministic (first pivot) | Randomized (2-way) | Randomized (3-way) |
|---:|---:|---:|---:|
| 1,000 | 14.0 ms (499,500) | 1.1 ms (11,042) | 1.3 ms (10,753) |
| 4,000 | 221.2 ms (7,998,000) | 5.1 ms (57,982) | 6.2 ms (54,604) |
| 16,000 | 3,730.4 ms (127,992,000) | 22.1 ms (257,048) | 28.6 ms (263,866) |
| 64,000 | skipped (> budget) | 100.8 ms (1,231,755) | 131.5 ms (1,212,189) |

**Reverse-sorted arrays**

| n | Deterministic (first pivot) | Randomized (2-way) | Randomized (3-way) |
|---:|---:|---:|---:|
| 1,000 | 25.1 ms (499,500) | 1.1 ms (11,014) | 1.3 ms (10,931) |
| 4,000 | 418.1 ms (7,998,000) | 5.1 ms (53,672) | 6.5 ms (55,111) |
| 16,000 | 6,532.4 ms (127,992,000) | 23.0 ms (265,758) | 29.1 ms (267,753) |
| 64,000 | skipped (> budget) | 100.2 ms (1,210,875) | 129.4 ms (1,201,591) |

**Repeated elements (only 10 distinct values)**

| n | Deterministic (first pivot) | Randomized (2-way) | Randomized (3-way) |
|---:|---:|---:|---:|
| 1,000 | 3.9 ms (53,646) | 4.3 ms (53,420) | 0.3 ms (3,487) |
| 4,000 | 57.5 ms (813,204) | 59.8 ms (814,843) | 1.0 ms (13,166) |
| 16,000 | 906.3 ms (12,864,071) | 930.0 ms (12,868,540) | 3.9 ms (49,858) |
| 64,000 | 14,495.1 ms (205,119,694) | 14,522.3 ms (205,015,130) | 17.8 ms (236,503) |

![Quicksort running times](results/quicksort_times.png)

![Comparisons on random input vs theory](results/quicksort_comparisons.png)

### 1.4 Discussion: relating the results to the theory

**1. Randomized quicksort shows n log n growth on every input with distinct keys.** When n grows by a factor of 4, an n log n algorithm's cost should grow by 4 · log(4n)/log(n). From n = 16,000 to 64,000 that factor is 4 × 15.97/13.97 ≈ 4.57. The measured factor for randomized 2-way on sorted input is 100.8/22.1 ≈ **4.56**. The random and reverse-sorted cases show the same behaviour. In the log-log plots the randomized curves are nearly identical for random, sorted, and reverse input. This matches the prediction that randomization removes the dependence on input order.

**2. Measured comparisons match the exact formula.** For n = 64,000, the formula 2(n+1)Hₙ − 4n predicts **1,234,438** comparisons. Randomized 2-way measured 1,218,612, which is 1.3% lower. For n = 8,000 the prediction is 121,051 and the measurement is 121,325. The number of comparisons in randomized quicksort has a standard deviation of about 0.65n, roughly 42,000 at n = 64,000. The observed differences are well within that range and are ordinary random variation. Normalised by n log₂ n, all curves move toward the limit 2 ln 2 ≈ 1.386, as predicted.

**3. Deterministic quicksort is exactly quadratic on sorted and reverse-sorted input.** It performed **exactly** n(n−1)/2 comparisons: for example, 127,992,000 = 16,000 × 15,999 / 2. Quadrupling n from 4,000 to 16,000 multiplied its time by 3,730/221 ≈ 16.9, close to the ×16 expected for Θ(n²). At n = 16,000 it was about **170× slower** than the randomized version on sorted input, and the gap keeps widening as n grows.

**4. Discrepancy: deterministic is slightly *faster* on random input.** Theory predicts the same expected number of comparisons for both, because a random array's first element is already a uniformly random pivot. The measured comparison counts do agree, within about 2%. Yet the deterministic version ran roughly 10–20% faster for n ≥ 4,000. The reason is overhead that the comparison model does not count. Each randomized partition calls `random.randint` and performs an extra swap, and in interpreted Python a call to `randint` is costly compared with a single loop iteration. This is a constant-factor difference, not an asymptotic one. It shows that the benefit of randomization is **protection against bad inputs**, not speed on inputs that are already random.

**5. Discrepancy: reverse-sorted is about 1.75–1.9× slower than sorted for deterministic quicksort, despite identical comparison counts.** Both perform exactly n(n−1)/2 comparisons. The difference lies in the swaps. On sorted input the pivot is the minimum, so the test `a[j] <= pivot` is never true and no swaps happen. On reverse-sorted input the first pivot is the maximum, so the test is true for every element and every element is swapped. The comparison count captures the asymptotic cost, but data movement changes the constant factor.

**6. Discrepancy: randomization does not help with heavy duplication when using a 2-way partition.** Both 2-way versions were quadratic on the repeated-elements input and had almost identical times. For example, both needed about 14.5 s at n = 64,000. This is not a contradiction of Section 1.2, because that proof assumes distinct keys. With 10 distinct values, the array contains about 10 blocks of n/10 equal keys. Each block costs about (n/10)²/2 comparisons, for a total of about n²/20. At n = 64,000 that gives **204.8 million**, and the measurement was **205.0 million**. The 3-way variant turns each block into a single pass. It ran about **800× faster** at n = 64,000, and its cost grew roughly in proportion to n (about 4× per quadrupling, with some run-to-run variation from the random pivots), as expected for O(n log k) with k fixed at 10.

**7. The price of 3-way partitioning.** On distinct keys the 3-way variant was about 20–30% slower than randomized 2-way. It can perform two key comparisons per element (`<`, then `>`), and its loop is more complex. Note that its counter records elements examined against the pivot, not individual `<`/`>` tests. This is the usual engineering trade-off: a small constant-factor cost on typical data in exchange for robustness on duplicate-heavy data. Production sorts make the same trade (Bentley & McIlroy, 1993).

**8. Measurement caveats.** Times below about 1 ms are affected by timer resolution and interpreter noise. The median of several runs reduces but does not remove that noise. Python's per-operation overhead also inflates every constant. A C implementation would be 50–100× faster in absolute terms, but its growth rates would be the same.

---

## Part 2: Hashing with Chaining

### 2.1 Implementation

`hash_table.py` implements `HashTableChaining`.

- **Hash function (universal family).** Each key is first mapped to an integer k: integers are used directly, and other keys go through Python's `hash()`. The value k is reduced modulo p = 2⁶¹ − 1, a Mersenne prime. The slot is then

    h_{a,b}(k) = ((a·k + b) mod p) mod m,  with a ∈ {1, …, p−1} and b ∈ {0, …, p−1} chosen at random.

  This is the Carter–Wegman family. For any two distinct keys k ≠ l in {0, …, p−1}, Pr_{a,b}[h(k) = h(l)] ≤ 1/m, so the family is **universal**. A new (a, b) pair is drawn every time the table is resized.
- **Chaining.** Each slot holds a Python list of `[key, value]` pairs.
- **Insert(key, value).** Hash the key, then scan its chain. If the key is already present, its value is replaced. Otherwise the pair is appended and the size is incremented. If the load factor then exceeds 0.75, the table doubles.
- **Search(key).** Hash the key, scan its chain, and return the value, or a default if the key is absent.
- **Delete(key).** Hash the key and scan its chain. When the key is found, it is overwritten with the chain's last pair and the last pair is popped. This removal takes O(1) because order within a chain does not matter. If the load factor then falls below 0.125, the table halves, but it never shrinks below its initial capacity.
- **Extras.** `len()`, `in`, `[]`, `load_factor`, `capacity`, `chain_lengths()`, and an option `resize=False` that turns off resizing for experiments.

`test_algorithms.py` checks the hash table against Python's built-in `dict` over 20,000 random insert, search, and delete operations, and it also checks resizing behaviour.

### 2.2 Expected running time under simple uniform hashing

**Assumption: simple uniform hashing (SUH).** Each key is equally likely to hash to any of the m slots, independently of all other keys. Let n be the number of stored keys and **α = n/m** the load factor. Computing h(k) takes O(1) time.

**Expected chain length.** Let nⱼ be the length of chain j. Then Σⱼ nⱼ = n, and by symmetry E[nⱼ] = n/m = α.

**Unsuccessful search.** A key that is not in the table hashes to a uniformly random slot, and the search must scan the entire chain. The expected number of keys examined is E[n_{h(k)}] = α, so the expected time is

  **Θ(1 + α).**

The "1" accounts for computing the hash and accessing the slot, which is necessary even when α is tiny.

**Successful search.** Suppose each of the n keys is equally likely to be the one searched for, and keys are inserted in the order x₁, …, xₙ. Because this implementation **appends** new keys to the end of their chain, the keys examined when searching for xᵢ are xᵢ itself plus every earlier key xⱼ (j < i) in the same slot. By SUH, Pr[h(xⱼ) = h(xᵢ)] = 1/m. The expected number examined is therefore

  (1/n) Σᵢ₌₁ⁿ ( 1 + Σⱼ₌₁ⁱ⁻¹ 1/m ) = 1 + (n−1)/(2m) = **1 + α/2 − α/(2n)**,

so the expected time is also **Θ(1 + α)**. (The CLRS proof inserts at the head of the chain and sums over later keys instead of earlier ones. The result is the same.)

**Insert.** Adding a pair to a chain is O(1). This implementation first scans the chain so that it never stores a duplicate key, so an insert costs the same as a search: expected **Θ(1 + α)**. Growth from resizing adds only **O(1) amortized** per insert (Section 2.4). An insert that allowed duplicates could skip the scan and run in O(1) worst case.

**Delete.** Deleting requires finding the key first, which is a successful search costing expected **Θ(1 + α)**. The removal itself is O(1) through the swap-and-pop technique. With a doubly linked list and a direct pointer to the node, deletion would be O(1) worst case. Shrinking adds O(1) amortized.

**Conclusion.** If the number of slots grows in proportion to the number of keys (m = Ω(n)), then α = O(1), and **insert, search, and delete all take O(1) expected time.**

**Universal hashing removes the SUH assumption.** SUH describes the keys, and real keys are rarely uniformly random. With a universal family, the randomness comes from the choice of hash function instead. For **any fixed set** of stored keys, the expected length of the chain that a key k maps to is at most α if k is absent and at most 1 + α if k is present (CLRS Theorem 11.3). The same Θ(1 + α) bounds therefore hold with no assumption about how the keys are distributed.

**Worst case.** If every key lands in one slot, every operation costs Θ(n). Universal hashing makes this extremely unlikely but not impossible. With α = 1 and a truly random hash function, the longest chain has expected length Θ(log n / log log n). The measured maximum chain length was 5 at α = 1 with n = 1,024.

### 2.3 How the load factor affects performance (empirical)

**Setup.** Resizing was disabled with m fixed at 1,024 slots. n = α·m random keys were inserted for each α. The search cost was then measured over 50,000 successful and 50,000 unsuccessful searches. Probe counts measure the number of keys compared in a chain, and times are the best of 5 runs.

| α | Successful: measured | Successful: theory 1+α/2−α/2n | Unsuccessful: measured | Unsuccessful: theory α | Longest chain | Hit time (µs) | Miss time (µs) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 1.13 | 1.12 | 0.25 | 0.25 | 3 | 0.38 | 0.36 |
| 0.5 | 1.29 | 1.25 | 0.49 | 0.50 | 4 | 0.40 | 0.38 |
| 0.75 | 1.36 | 1.37 | 0.76 | 0.75 | 4 | 0.40 | 0.39 |
| 1 | 1.50 | 1.50 | 1.01 | 1.00 | 5 | 0.40 | 0.40 |
| 2 | 1.98 | 2.00 | 1.99 | 2.00 | 7 | 0.43 | 0.43 |
| 4 | 2.99 | 3.00 | 4.02 | 4.00 | 12 | 0.47 | 0.50 |
| 8 | 5.06 | 5.00 | 7.99 | 8.00 | 19 | 0.52 | 0.59 |
| 16 | 9.01 | 9.00 | 16.04 | 16.00 | 31 | 0.74 | 1.04 |

![Search cost vs load factor](results/hash_load_factor.png)

**Observations**

- The measured probe counts match the formulas within about 3% at every load factor. Unsuccessful searches cost about α probes, and successful searches about 1 + α/2. Both grow **linearly in α**, as the Θ(1 + α) bound predicts.
- For α < 2, successful searches cost *more* than unsuccessful ones. A successful search always examines at least the key it finds, while an unsuccessful search in an empty slot examines nothing. Above α = 2, unsuccessful searches cost more, because they must scan the whole chain while a successful search stops, on average, halfway.
- Wall-clock time rises more slowly than the probe count. Each search has a fixed overhead for the function call, the hash computation, and the Python loop setup. That overhead corresponds to the "1" in Θ(1 + α) and dominates at small α. At α = 16, the probe count is 64× higher than at α = 0.25, but time is only about 2.9× higher. At large enough α the α term dominates, as the fixed-table results in Section 2.5 show.
- The longest chain grows with α as well. At α = 16 some searches scanned 31 keys, nearly twice the average. A high load factor therefore also increases the variance of operation times.

### 2.4 Keeping the load factor low and minimizing collisions

**1. Dynamic resizing (implemented).**

- **Growth rule.** When an insert would push α above α_max = 0.75, allocate a table with 2m slots, draw a new hash function, and reinsert every key. Rehashing costs Θ(n), but it happens only after the number of keys has doubled since the previous rehash.
- **Amortized cost of growth.** Starting from capacity m₀ and growing to hold n keys, the rehashes move at most α_max·(m₀ + 2m₀ + 4m₀ + … + m_final) < 2·α_max·m_final = O(n) keys in total. Spread over n inserts, that is **O(1) amortized** per insert.
- **Shrink rule.** When a delete pushes α below α_min = 0.125, halve the table. This keeps memory proportional to n after mass deletions.
- **Hysteresis.** α_min is set well below α_max/2. Right after growth α ≈ 0.375, and right after shrinking α ≈ 0.25, so Ω(n) operations must occur before another resize. If the two thresholds were adjacent, for example growing at 0.75 and shrinking at 0.375, an alternating insert/delete sequence at the boundary would trigger a Θ(n) rehash on every operation. With the gap, deletions are also O(1) amortized.
- **Evidence.** In Experiment B the resizing table ended every run at α ≈ 0.49, regardless of n.

**2. A good hash function (implemented).**

- A universal family provides the 1/m collision bound for **every** key set, including structured, non-random key sets that break simple functions.
- Experiment C illustrates this with 4,096 keys that are all multiples of m = 1,024:

| Hash function | Colliding pairs | Longest chain | Slots used |
|---|---:|---:|---:|
| Division method h(k) = k mod m | 8,386,560 (every pair) | 4,096 | 1 |
| Universal family, averaged over 100 random draws | 7,497 (bound C(n,2)/m = 8,190) | median 5, worst 24 | median 1,024 |

- The division method puts every key into slot 0 and turns the table into a linked list. The universal family averages below its theoretical bound of C(n,2)/m colliding pairs.
- The worst single draw still produced a chain of 24. Universality is a guarantee **on average over the random choice of function**, not for each individual function. Drawing a new function on every resize ensures that an unlucky draw is not permanent.

**3. Other strategies (not implemented, but relevant):**

- **Pre-sizing.** If the final n is known in advance, choosing m ≈ n/α_max avoids all intermediate rehashes.
- **Incremental resizing.** A single Θ(n) rehash causes a latency spike. Real-time systems avoid this by keeping both the old and new tables and moving a few entries during each operation, as Redis does.
- **Bounding chain length.** Java's `HashMap` converts a chain into a balanced binary search tree once the chain grows past 8 entries. This caps the worst case at O(log n) even when many keys collide.
- **Avoiding weak moduli.** With the plain division method, m should be a prime not close to a power of 2. Otherwise m = 2ʳ selects only the low bits of the key. The multiply-mod-p step in the universal family avoids this problem.
- **Choosing α_max.** A lower threshold, such as 0.5, gives shorter chains at the cost of more memory. A higher one, such as 1–2, saves memory at the cost of slower operations. For chaining, values between about 0.75 and 1 are typical.

### 2.5 Scalability with and without resizing (empirical)

**Setup.** Experiment B inserts n random keys, searches for all of them in shuffled order, and then deletes them all. The table reports time per operation in microseconds.

| Mode | n | Final capacity | Final α | Insert (µs) | Search (µs) | Delete (µs) |
|---|---:|---:|---:|---:|---:|---:|
| resizing | 1,000 | 2,048 | 0.49 | 1.51 | 0.40 | 1.14 |
| resizing | 4,000 | 8,192 | 0.49 | 2.20 | 0.42 | 1.20 |
| resizing | 16,000 | 32,768 | 0.49 | 1.84 | 0.53 | 1.38 |
| resizing | 64,000 | 131,072 | 0.49 | 2.58 | 0.65 | 1.78 |
| resizing | 256,000 | 524,288 | 0.49 | 3.82 | 1.13 | 2.92 |
| resizing | 1,000,000 | 2,097,152 | 0.48 | 6.32 | 1.26 | 4.36 |
| fixed m = 1,024 | 1,000 | 1,024 | 0.98 | 0.59 | 0.53 | 0.63 |
| fixed m = 1,024 | 4,000 | 1,024 | 3.91 | 0.56 | 0.46 | 0.66 |
| fixed m = 1,024 | 16,000 | 1,024 | 15.62 | 0.91 | 0.68 | 0.81 |
| fixed m = 1,024 | 64,000 | 1,024 | 62.50 | 2.14 | 2.19 | 1.73 |
| fixed m = 1,024 | 256,000 | 1,024 | 250.00 | 8.66 | 9.36 | 6.08 |

![Hash table scaling](results/hash_scaling.png)

**Discussion**

- **The fixed table degrades linearly.** As n grows 256× from 1,000 to 256,000, α grows 256× as well, and search time per operation rises from 0.53 µs to 9.36 µs, about 18×. The growth is slower than 256× because of the fixed per-operation overhead described in Section 2.3. From n = 64,000 to 256,000, however, α quadruples and search time grows about 4.3×, so the α term clearly dominates there. Over a whole workload, a table without resizing costs Θ(n²/m) = Θ(n²) in total.
- **The resizing table keeps the number of probes constant.** The load factor stays at about 0.5 for every n, so the expected chain length does not change. Search time grew only from 0.40 µs to 1.26 µs over a **1,000× increase in n**, compared with the 18× increase for the fixed table over just 256× growth.
- **Discrepancy: why the resizing table's time is not perfectly flat.** The RAM model used in Section 2.2 assumes every memory access costs the same. On real hardware that is not true. At n = 1,000,000 the table holds about 2 million Python list objects, many megabytes spread across the heap. Most accesses then miss the CPU cache, while at n = 1,000 the entire table fits in cache. The probe count stays constant, but each probe becomes more expensive. This effect is a property of the memory hierarchy, not of the algorithm. It appears in any large data structure, including the quicksort arrays in Part 1.
- **Insert and delete cost more than search in the resizing table.** Their per-operation times include the amortized rehashing work: growth for inserts and shrinking for deletes. Each rehash also allocates new list objects in Python. The fixed table is cheaper at very small n because it never rehashes. This is the cost that dynamic resizing trades for scalability, and the trade clearly pays off as soon as α would otherwise exceed a small constant.

---

## Conclusions and algorithm selection

1. **Randomized quicksort** gives expected Θ(n log n) time on every input of distinct keys. The measured comparisons matched 2(n+1)Hₙ − 4n within ordinary random variation. On random data it is roughly 10–20% slower than deterministic quicksort because of the cost of generating random numbers. On sorted or reverse-sorted data it is over 100× faster at n = 16,000, and the gap grows with n. When the input order is unknown or could be adversarial, the pivot should be randomized.
2. **Pivot randomization alone does not handle duplicates.** With a 2-way partition, heavily duplicated data is quadratic whether or not the pivot is random. For real data, which often contains repeated keys, **randomized quicksort with 3-way partitioning** is the robust choice. It was about 800× faster on duplicate-heavy input and only 20–30% slower on distinct keys.
3. **Hashing with chaining** supports insert, search, and delete in Θ(1 + α) expected time. The experiments confirmed the exact constants α and 1 + α/2. Keeping α bounded requires dynamic resizing by doubling, which adds only O(1) amortized cost per operation. A universal hash function, re-drawn on every resize, keeps these guarantees independent of the key distribution.
4. **General lesson.** Asymptotic analysis correctly predicted every growth rate observed in these experiments. The constant factors came from effects outside the model, such as interpreter overhead, random number generation, data movement, and cache behaviour. Both kinds of analysis are needed to choose an algorithm.

## References

- Bentley, J. L., & McIlroy, M. D. (1993). Engineering a sort function. *Software: Practice and Experience, 23*(11), 1249–1265.
- Carter, J. L., & Wegman, M. N. (1979). Universal classes of hash functions. *Journal of Computer and System Sciences, 18*(2), 143–154.
- Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2022). *Introduction to algorithms* (4th ed.). MIT Press. (Chapter 7: Quicksort; Chapter 11: Hash Tables; Chapter 16: Amortized Analysis.)
- Hoare, C. A. R. (1962). Quicksort. *The Computer Journal, 5*(1), 10–16.
