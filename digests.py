from collections import Counter

tolerance = 0
Digest = tuple[list[str], list[float]]

class f:
    def __init__(self, start = None, end = None, length = 0, position = 0):
        self.startEnzyme = start #there is redundancy here, but I like it
        #this means not needing to find the enzyme responsible for the start of the fragment
        self.endEnzyme = end #enzyme responsible for end of the fragment
        self.length = length #length of fragment
        self.startPosition = position #start position of fragment

#n_sum, the first helper function, mainly replaced by coverage_possibilities
def n_sum(fragments, total):
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
    sampling = []
    for digest in chain:
        sampling.append(sum(digest[1]))
    average = sum(sampling)/len(sampling)
    for length in sampling:
        if not (0.9*average < length and length < 1.1*average):
            return False
    return average

def build_canonical(chain):
    if not find_length(chain):
        print("Bad digests")
        return 
    plasmid = [f(None, None, find_length(chain))]
    ndigest = recurse(plasmid, chain, chain)
    return ndigest

def simulate_digest(current_map, enzymes):
    """
    Simulates a restriction digest on the current map using specified enzymes.
    Returns a sorted list of fragment lengths.
    """
    #collect all cut positions for these enzymes
    cut_positions = set()
    
    for fragment in current_map:
        #check if this fragment's end is cut by one of our enzymes
        if fragment.endEnzyme in enzymes:
            end_position = fragment.startPosition + fragment.length
            cut_positions.add(end_position)
    
    #convert to sorted list
    cut_positions = sorted(cut_positions)
    
    if not cut_positions:
        # No cuts, return total length
        return [sum(frag.length for frag in current_map)]
    
    #for circular plasmid, calculate fragments between cuts
    total_length = sum(frag.length for frag in current_map)
    simulated_frags = []
    
    for i in range(len(cut_positions)):
        start = cut_positions[i]
        end = cut_positions[(i + 1) % len(cut_positions)]
        
        if i == len(cut_positions) - 1:
            #last fragment wraps around
            simulated_frags.append(total_length - start + cut_positions[0])
        else:
            simulated_frags.append(end - start)
    
    return sorted(simulated_frags)


def check_validity(current_map, chain):
    """
    Validates the current mapping against all digests in the chain.
    """
    for digest in chain:
        enzymes, expected_frags = digest
        simulated = simulate_digest(current_map, enzymes)
        
        if Counter(simulated) != Counter(expected_frags):
            return False
    
    return True


def find_mismatches(simulated_frags, actual_frags):
    """
    Finds fragments that don't match between simulated and actual digests.
    Returns (mismatched_old, mismatched_new) where old fragments were cut into new fragments.
    """
    sim_counter = Counter(simulated_frags)
    actual_counter = Counter(actual_frags)
    
    #find fragments only in simulated (these got cut)
    mismatched_old = []
    for frag, count in sim_counter.items():
        actual_count = actual_counter.get(frag, 0)
        if count > actual_count:
            mismatched_old.extend([frag] * (count - actual_count))
    
    #find fragments only in actual (these are the cut pieces)
    mismatched_new = []
    for frag, count in actual_counter.items():
        sim_count = sim_counter.get(frag, 0)
        if count > sim_count:
            mismatched_new.extend([frag] * (count - sim_count))
    
    return sorted(mismatched_old, reverse=True), sorted(mismatched_new, reverse=True)


def insert_cuts(current_map, division, old_frags, new_enzyme, known_enzymes):
    """
    Inserts new enzyme cuts into the map based on the division pattern.
    division: list of lists, where each sublist shows how an old fragment was divided
    old_frags: the fragments that were cut
    new_enzyme: the enzyme making the new cuts
    
    Example: A 10000 A -> A 4000 B | B 3000 B | B 2000 B | B 1000 A
    """
    new_map = []
    old_frag_iter = iter(old_frags)
    current_old_frag = next(old_frag_iter, None)
    division_index = 0
    
    for fragment in current_map:
        #checks if this fragment needs to be divided
        if current_old_frag is not None and abs(fragment.length - current_old_frag) < tolerance + 0.01:
            #gets divided according to division pattern
            sub_frags = division[division_index]
            position = fragment.startPosition
            
            for i, sub_length in enumerate(sub_frags):
                if i == 0:
                    #first piece: keeps original start enzyme, ends with new enzyme
                    new_map.append(f(fragment.startEnzyme, new_enzyme, sub_length, position))
                elif i < len(sub_frags) - 1:
                    #middle pieces: start and end with new enzyme
                    new_map.append(f(new_enzyme, new_enzyme, sub_length, position))
                else:
                    #last piece: starts with new enzyme, keeps original end enzyme
                    new_map.append(f(new_enzyme, fragment.endEnzyme, sub_length, position))
                position += sub_length
            
            division_index += 1
            current_old_frag = next(old_frag_iter, None)
        else:
            #fragment unchanged, just maybe update position
            new_map.append(f(fragment.startEnzyme, fragment.endEnzyme, fragment.length, fragment.startPosition))
    return new_map


def recurse(current_map, remaining_digests, chain):
    """
    Recursively builds the plasmid map by placing enzyme cut sites.
    Compares simulated digest with actual digest to find where new enzymes cut.
    """
    visualize_map(current_map)
    
    if not remaining_digests:
        return current_map if check_validity(current_map, chain) else None
    
    current_digest = remaining_digests[0]
    new_enzymes, new_fragments = current_digest
    
    known_enzymes = set()
    for frag in current_map:
        if frag.endEnzyme is not None:
            known_enzymes.add(frag.endEnzyme)
        if frag.startEnzyme is not None:
            known_enzymes.add(frag.startEnzyme)
    
    new_enzyme_set = set(new_enzymes) - known_enzymes
    
    #first case: First digest on uncut plasmid
    if not known_enzymes and len(current_map) == 1:
        print(f"First digest: {new_enzymes}, fragments: {new_fragments}")
        
        #for a single enzyme on circular plasmid, all fragments connect via that enzyme
        #try all rotations of the fragments
        for rotation in range(len(new_fragments)):
            rotated = new_fragments[rotation:] + new_fragments[:rotation]
            
            new_map = []
            position = 0
            enzyme = new_enzymes[0]
            
            for i, length in enumerate(rotated):
                #all cuts by the same enzyme in circular plasmid
                new_map.append(f(enzyme, enzyme, length, position))
                position += length
            
            if check_validity(new_map, chain[:1]):
                result = recurse(new_map, remaining_digests[1:], chain)
                if result is not None:
                    return result
        
        return None
    
    if not new_enzyme_set:
        if check_validity(current_map, [current_digest]):
            return recurse(current_map, remaining_digests[1:], chain)
        return None
    
    new_enzyme = list(new_enzyme_set)[0]
    simulated_frags = simulate_digest(current_map, known_enzymes)
    mismatched_old, mismatched_new = find_mismatches(simulated_frags, new_fragments)
    
    if not mismatched_old:
        if check_validity(current_map, chain[:chain.index(current_digest) + 1]):
            return recurse(current_map, remaining_digests[1:], chain)
        return None
    
    possible_divisions = coverage_possibilities(mismatched_new, mismatched_old)
    
    for division in possible_divisions:
        new_map = insert_cuts(current_map, division, mismatched_old, new_enzyme, known_enzymes)
        
        if new_map and check_validity(new_map, chain[:chain.index(current_digest) + 1]):
            result = recurse(new_map, remaining_digests[1:], chain)
            if result is not None:
                return result
    
    return None

def visualize_map(current_map):
    total_length = sum(frag.length for frag in current_map)
    print(f"Total length: {total_length}; Number of fragments: {len(current_map)}")
    print()
    
    for i, frag in enumerate(current_map):
        start_label = frag.startEnzyme if frag.startEnzyme else "START"
        end_label = frag.endEnzyme if frag.endEnzyme else "END"
        print(f"  Fragment {i}: {start_label} --[{frag.length}bp]--> {end_label}  (pos: {frag.startPosition})")
    repr_str = ""
    #show linear representation
    for frag in current_map:
        start = frag.startEnzyme if frag.startEnzyme else "**"
        repr_str += f"{start}={frag.length}="
    end = current_map[-1].endEnzyme if current_map[-1].endEnzyme else "**"
    repr_str += end
    print(f"  {repr_str}")
    print(f"{'='*60}\n")

A = (["EcoRI"], [6000, 4000])
B = (["BamHI"], [5000, 3000, 2000])
C = (["HindIII"], [7000, 3000])
D = (["EcoRI", "BamHI"], [3000,2000,2000,1000,2000])
E = (["BamHI", "HindIII", "PstI"], [3000,2000,2000,1500,1500])
F = (["EcoRI", "HindIII", "XhoI"], [3000,2500,2000,1500,1000])
#test = [A, B, C, D, E, F]
test = [A, B, D]

chain = build_chain(test)
print(chain)
mapping = build_canonical(chain)
if mapping:
    for frag in mapping:
        print(frag)
else:
    print("Something either went wrong OR your digests are invalid.")

#print(coverage_possibilities(D[1], A[1]))
#print(coverage_possibilities(D[1], B[1]))


