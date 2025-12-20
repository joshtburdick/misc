
import statistics
import random

def solve():
    # Fuzzing to verify m_B > m_A implies strict inequality
    for _ in range(1000):
        n = random.randint(1, 20)
        k = random.randint(1, 50)
        
        # Create random distinct A
        # To avoid easy cases, make A tight or spread
        A = sorted(random.sample(range(0, 1000), n))
        m_A = statistics.median(A)
        
        # Construct minimal B
        # B must contain k distinct values > a_i for each i
        # To minimize m_B, we want B values as small as possible.
        # However, they must match distinctly.
        # We can greedy assign from smallest a_i? No, smallest a_i allows smallest b.
        # But we need disjoint sets.
        # Constraints:
        # S_1 > a_1
        # ...
        # S_n > a_n
        # Optimal strategy to minimize B values:
        # Just pick them?
        # Let's just generate a VALID B randomly and check. 
        # Actually random check is weaker than "minimal" check.
        # But if the theorem holds for "all" B, it holds for minimal B.
        
        # Let's generate a minimal valid B.
        # We need kn values.
        # We need k values > a_n
        # We need k values > a_{n-1} (distinct from above)
        # ...
        # The tightest packing:
        # assigned_B = set()
        # processing from largest a_i down to a_1 allows us to pick largest needed values first?
        # No, to MINIMIZE median, we want values as small as possible.
        # So we should satisfy a_1 (smallest) with smallest possible values?
        # But those values must be > a_1.
        # Wait, values for a_n MUST be > a_n.
        # Values for a_1 only need be > a_1.
        # So we should use larger numbers for harder constraints (a_n).
        # Actually, let's just pick k values for a_1: a_1+1 ... a_1+k (checking collisions)
        # Then k for a_2: a_2+1 ... (checking collisions)
        # This greedy strategy works?
        
        used_B = set()
        B_list = []
        for x in A:
            # Find k smallest integers > x not in used_B
            found_count = 0
            curr = x + 1
            while found_count < k:
                if curr not in used_B:
                    used_B.add(curr)
                    B_list.append(curr)
                    found_count += 1
                curr += 1
        # print(A)
        # print(sorted(B_list))
        m_B = statistics.median(B_list)
        diff = statistics.median(B_list) - statistics.median(A)
        print(f"m_A: {m_A}, m_B: {m_B}, k: {k},diff: {diff}")
        if m_B <= m_A:
            print(f"COUNTEREXAMPLE FOUND!")
            print(f"A: {A}")
            print(f"B: {sorted(B_list)}")
            print(f"m_A: {m_A}, m_B: {m_B}")
            return

    print("No counterexamples found in 1000 iterations.")

solve()
