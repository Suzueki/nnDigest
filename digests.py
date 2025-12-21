tolerance = 0

Digest = tuple[list[str], list[float]]

class f:
    def __init__(self, start = None, end, length):
        self.start = start
        self.end = end
        self.length = length

#n-sum, should help in the majority of digests
def nsum(fragments, total):
    groups = []
    n = len(fragments)

    def backtrack(x, remain, current_path):
        if abs(remain) <= tolerance:
            groups.append(current_path.copy())
            return

        if remain < -tolerance:
            return

        for i in range (x, n):
            if remain - fragments[i] >= -tolerance:
                current_path.append(i)
                #becomes an n-1-sum problem to find the others in the group.
                #much better than the 2^n memory or n^2 time approach
                backtrack(i+1, remain-fragments[i], current_path)
                current_path.pop()
    backtrack(0, total, [])
    return groups

def coverage_possibilities(subfragments, superfragments):
    subfragments = sorted(subfragments, reverse = True)
    superfragments = sorted(superfragments, reverse = True)

    results = []
    def backtrack(remaining_frags, superfrag_index, curr):
        if superfrag_index == len(superfragments):
            if not remaining_frags:
                results.append([list(f) for f in curr])
            return

        capacity = superfragments[superfrag_index]

        def choose_frags(start, curr_frags, curr_length):
            if curr_length == capacity:
                yield list(curr_frags)
                return
            if curr_length > capacity:
                return
            prev = None
            for i in range(start, len(remaining_frags)):
                val = remaining_frags[i]
                if val == prev or curr_length + val > capacity:
                    continue
                prev = val
                curr_frags.append(val)
                yield from choose_frags(i + 1, curr_frags, curr_length + val)
                curr_frags.pop()
                
        used_subsets = set()        
        for subset in choose_frags(0, [], 0):
            subset_key = tuple(sorted(subset))
            if subset_key in used_subsets:
                continue
            used_subsets.add(subset_key)
            remaining = remaining_frags.copy()
            for x in subset:
                remaining.remove(x)
            backtrack(
                remaining, superfrag_index + 1, curr + [subset]
            )
    backtrack(subfragments, 0, [])
    return results
                
def find_boundaries(digest1, digest2):
    #quick for perfect accuracy, should be good for manual inputs
    if tolerance == 0:
        return set(digest1) & set(digest2)

    #log() approach, should work if bad manual inputs or read from a picture
    shared = set()
    sorted1 = sorted(digest1)
    sorted2 = sorted(digest2)
    shared = set()
    i = 0
    j = 0

    #2 pointer approach
    while i < len(sorted1) and j < len(sorted2):
        diff = sorted1[i] - sorted2[j]
        if abs(diff) <= tolerance:
            shared.add(sorted1[i])
            i = i + 1
            j = j + 1
        elif diff < 0:
            i = i + 1
        else:
            j = j + 1
    return shared


def build_superfragments(superfrags, subfrags):
    groups = []
    for superfrag in superfrags:
        subgroups = nsum(subfrags, superfrag)
        groups.append((superfrag, subgroups)) 
    return groups

def build_chain(digests):
    seen_enzymes = set()
    seen_fragments = set()
    chain = []
    remaining = digests.copy()

    while remaining:
        scored = []
        for d in remaining:
            enzymes, frags = d
            new_enzymes = set(enzymes) - seen_enzymes
            new_frags = set(frags) - seen_fragments
            common_enzymes = set(enzymes) & seen_enzymes
            #score = (len(new_enzymes), len(new_frags))
            score = (len(new_enzymes) + 1)**2 + len(new_frags) - ((len(common_enzymes) + 3)**1.5)
            scored.append((score, d))
        #choose digest with minimal (new enzymes, new fragments)
        scored.sort(key=lambda x: x[0])
        best = scored[0][1]
        #update state
        enzymes, frags = best
        seen_enzymes.update(enzymes)
        seen_fragments.update(frags)
        chain.append(best)
        remaining.remove(best)
    return chain

def find_length(chain):
    

def build_canonical(chain):
    plasmid = f(None, None, find_length(chain))

A = (["EcoRI"], [6000, 4000])
B = (["BamHI"], [5000, 3000, 2000])
C = (["HindIII"], [7000, 3000])
D = (["EcoRI", "BamHI"], [3000,2000,2000,1000,2000])
E = (["BamHI", "HindIII", "PstI"], [3000,2000,2000,1500,1500])
F = (["EcoRI", "HindIII", "XhoI"], [3000,2500,2000,1500,1000])
test = [A, B, C, D, E, F]

print(build_chain(test))

#print(coverage_possibilities(D[1], A[1]))
#print(coverage_possibilities(D[1], B[1]))


