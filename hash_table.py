"""
Hash table with separate chaining and universal hashing (Assignment 3, Part 2).

Hash function
-------------
Keys are first turned into a non-negative integer k (integers are used as-is,
other hashable keys go through Python's built-in hash()). The slot is then
chosen with the Carter-Wegman universal family

        h_{a,b}(k) = ((a * k + b) mod p) mod m

where p = 2^61 - 1 is a Mersenne prime, a is drawn uniformly from {1..p-1},
b uniformly from {0..p-1}, and m is the number of slots. For any two distinct
keys k != l (with k, l < p) the probability over the random choice of (a, b)
that they collide is at most 1/m. A new (a, b) pair is drawn every time the
table is resized.

Collision resolution
--------------------
Each slot holds a Python list of [key, value] pairs (the "chain").

Dynamic resizing
----------------
  * grow  : capacity doubles when the load factor would exceed MAX_LOAD (0.75)
  * shrink: capacity halves when the load factor drops below MIN_LOAD (0.125),
            never going below the initial capacity.
Doubling / halving keeps the amortized cost of insert and delete at O(1).
"""

import random

_MERSENNE_PRIME = (1 << 61) - 1


class HashTableChaining:
    MAX_LOAD = 0.75
    MIN_LOAD = 0.125

    def __init__(self, capacity=8, resize=True, seed=None):
        if capacity < 1:
            raise ValueError("capacity must be at least 1")
        self._rng = random.Random(seed)
        self._initial_capacity = capacity
        self._resize_enabled = resize
        self._size = 0
        self._m = capacity
        self._table = [[] for _ in range(capacity)]
        self._choose_hash_function()

    # ------------------------------------------------------------------ hashing
    def _choose_hash_function(self):
        """Draw a fresh function from the universal family."""
        p = _MERSENNE_PRIME
        self._a = self._rng.randrange(1, p)
        self._b = self._rng.randrange(0, p)

    def _hash(self, key):
        k = key if isinstance(key, int) and not isinstance(key, bool) else hash(key)
        k %= _MERSENNE_PRIME  # map into the universe {0..p-1}
        return ((self._a * k + self._b) % _MERSENNE_PRIME) % self._m

    # ------------------------------------------------------------ core methods
    def insert(self, key, value):
        """Add key -> value. If key already exists its value is replaced."""
        chain = self._table[self._hash(key)]
        for pair in chain:
            if pair[0] == key:
                pair[1] = value
                return
        chain.append([key, value])
        self._size += 1
        if self._resize_enabled and self._size / self._m > self.MAX_LOAD:
            self._rehash(self._m * 2)

    def search(self, key, default=None):
        """Return the value stored for key, or `default` if key is absent."""
        for k, v in self._table[self._hash(key)]:
            if k == key:
                return v
        return default

    def delete(self, key):
        """Remove key. Returns True if it was present, False otherwise."""
        chain = self._table[self._hash(key)]
        for i, pair in enumerate(chain):
            if pair[0] == key:
                # Order inside a chain does not matter: swap with last, pop (O(1)).
                chain[i] = chain[-1]
                chain.pop()
                self._size -= 1
                if (self._resize_enabled
                        and self._m > self._initial_capacity
                        and self._size / self._m < self.MIN_LOAD):
                    self._rehash(max(self._initial_capacity, self._m // 2))
                return True
        return False

    # ----------------------------------------------------------------- resizing
    def _rehash(self, new_capacity):
        old_items = [pair for chain in self._table for pair in chain]
        self._m = new_capacity
        self._table = [[] for _ in range(new_capacity)]
        self._choose_hash_function()
        for pair in old_items:
            self._table[self._hash(pair[0])].append(pair)

    # ------------------------------------------------------------ introspection
    @property
    def capacity(self):
        return self._m

    @property
    def load_factor(self):
        return self._size / self._m

    def chain_lengths(self):
        return [len(c) for c in self._table]

    def __len__(self):
        return self._size

    def __contains__(self, key):
        sentinel = object()
        return self.search(key, sentinel) is not sentinel

    def __getitem__(self, key):
        sentinel = object()
        v = self.search(key, sentinel)
        if v is sentinel:
            raise KeyError(key)
        return v

    def __setitem__(self, key, value):
        self.insert(key, value)

    def __delitem__(self, key):
        if not self.delete(key):
            raise KeyError(key)

    def items(self):
        for chain in self._table:
            for k, v in chain:
                yield k, v

    def __repr__(self):
        return (f"HashTableChaining(size={self._size}, capacity={self._m}, "
                f"load_factor={self.load_factor:.3f})")


if __name__ == "__main__":
    t = HashTableChaining(seed=42)
    for word in ["apple", "banana", "cherry", "date", "elderberry"]:
        t.insert(word, len(word))
    print(t)
    print("search('cherry') ->", t.search("cherry"))
    print("delete('banana') ->", t.delete("banana"))
    print("search('banana') ->", t.search("banana"))
    print(t)
