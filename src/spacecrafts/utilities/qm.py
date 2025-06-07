"""Implementation of the Quine McCluskey algorithm for boolean minimization."""

from typing import Iterable, List, Set, Tuple


def _combine_terms(term1: str, term2: str) -> Tuple[str, bool]:
    """Combine two terms if they differ by one bit.

    Returns (combined_term, combined?)"""
    diff = None
    for i, (a, b) in enumerate(zip(term1, term2)):
        if a != b:
            if diff is not None:
                return term1, False
            diff = i
    if diff is None:
        return term1, False
    return term1[:diff] + '-' + term1[diff + 1:], True


def _initial_group(terms: Set[int], num_vars: int):
    groups = {}
    for t in terms:
        b = format(t, f"0{num_vars}b")
        groups.setdefault(b.count('1'), set()).add(b)
    return groups


def _iterate_combine(groups: dict, num_vars: int):
    new_groups = {}
    checked = set()
    prime = set()
    keys = sorted(groups.keys())
    for i, k in enumerate(keys[:-1]):
        g1 = groups[k]
        g2 = groups[keys[i + 1]]
        for term1 in g1:
            combined_flag = False
            for term2 in g2:
                combined, ok = _combine_terms(term1, term2)
                if ok:
                    new_groups.setdefault(combined.count('1'), set()).add(combined)
                    combined_flag = True
            if not combined_flag:
                prime.add(term1)
        # terms in last group may become prime later
    for k in keys[-1:]:
        for term in groups[k]:
            if all(_combine_terms(term, other)[1] is False for other in groups.get(k + 1, [])):
                prime.add(term)
    return new_groups, prime


def _get_prime_implicants(terms: Set[int], dontcares: Set[int], num_vars: int) -> Set[str]:
    all_terms = terms | dontcares
    groups = _initial_group(all_terms, num_vars)
    prime_implicants = set()
    while groups:
        groups, primes = _iterate_combine(groups, num_vars)
        prime_implicants.update(primes)
    return prime_implicants


def _covers(term: str, number: int) -> bool:
    bits = format(number, f"0{len(term)}b")
    return all(t == '-' or t == b for t, b in zip(term, bits))


def _essential_primes(prime_implicants: Set[str], minterms: Set[int]) -> Tuple[Set[str], Set[int]]:
    chart = {m: set() for m in minterms}
    for m in minterms:
        for p in prime_implicants:
            if _covers(p, m):
                chart[m].add(p)
    essential = set()
    for m, implicants in chart.items():
        if len(implicants) == 1:
            essential.update(implicants)
    remaining_minterms = {m for m in minterms if not any(_covers(p, m) for p in essential)}
    return essential, remaining_minterms


def _petrick_method(prime_implicants: Set[str], minterms: Set[int]) -> Set[str]:
    chart = {m: {p for p in prime_implicants if _covers(p, m)} for m in minterms}
    P = [ {frozenset([p]) for p in ps} for ps in chart.values() ]
    while len(P) > 1:
        A = P.pop()
        B = P.pop()
        new = set()
        for a in A:
            for b in B:
                new.add(a | b)
        # remove supersets
        minimal = set()
        for c in new:
            if not any(c > d for d in new if c != d):
                minimal.add(c)
        P.append(minimal)
    solutions = P[0]
    best = min(solutions, key=lambda x: sum(len(p.replace('-', '')) for p in x))
    return set(best)


def quine_mccluskey(minterms: Iterable[int], dontcares: Iterable[int] = None, num_vars: int = None) -> List[str]:
    """Minimize boolean function using Quine McCluskey algorithm.

    Parameters
    ----------
    minterms : iterable of int
        Minterms that should evaluate to 1.
    dontcares : iterable of int, optional
        Terms that can be ignored in the minimization.
    num_vars : int, optional
        Number of variables. If not provided it is inferred from the
        maximum term.

    Returns
    -------
    list of str
        Simplified implicants represented with '-' as don't care.
    """
    minterms = set(minterms)
    dontcares = set(dontcares or [])
    if not minterms and not dontcares:
        return []
    if num_vars is None:
        num_vars = max(minterms | dontcares).bit_length()
    primes = _get_prime_implicants(minterms, dontcares, num_vars)
    essential, remaining = _essential_primes(primes, minterms)
    if remaining:
        cover = _petrick_method(primes - essential, remaining)
    else:
        cover = set()
    result = sorted(essential | cover)
    return result
