import copy
from collections import Counter
from itertools import product

tol = 0
all_solutions = False
solutions = []

'''
Infinite randomized test cases implemented tomorrow.
'''

class Mapping:
    def __init__(self, length = 0):
        self.sites = {}
        self.length = length
        self.known_enzymes = set()

    def add_site(self, enzymes, position):
        #now handles both a single string "EcoRI" or a set {"EcoRI", "BamHI"}!
        #extremely important for sites cleaved by multiple enzymes in an n-digest
        if isinstance(enzymes, str):
            enzymes = {enzymes}
        pos = round(float(position) % self.length, 4)
        if pos not in self.sites:
            self.sites[pos] = set()
        self.sites[pos].update(enzymes)
        self.known_enzymes.update(enzymes)

    def add_sites(self, cuts):
        #uses add_site extensively
        #cuts is an list of (position, enzyme) tuples
        #we always permute outside the loop because it's essential for recursion
        for position, enzyme in cuts:
            self.add_site(enzyme, position)

    def get_buckets(self, active_enzymes):
        """
        Finds the 'buckets/contarios' for the current digest.
        Returns a list of (start_pos, length) tuples.
        """
        #get positions of sites that contain any of the enzymes in the current digest
        active_pos = self.get_sorted_positions(active_enzymes)
        #if none of these enzymes are on the map yet (seed map at pos 0)
        if not active_pos:
            return [(0.0, float(self.length))]
        buckets = []
        for i in range(len(active_pos)):
            p1 = active_pos[i]
            p2 = active_pos[(i + 1) % len(active_pos)]
            #allows for circular plasmids, essentially a modulo operation
            dist = (p2 - p1) if p2 > p1 else (self.length - p1 + p2)
            buckets.append((p1, round(dist, 4)))
        return buckets

    def get_active_sites_sorted(self, enzymes):
        """Helper to get (pos, enzymes) tuples for specific enzymes, sorted by pos."""
        active_pos = [(pos, enz_set) for pos, enz_set in self.sites.items() if not enz_set.isdisjoint(enzymes)]
        return sorted(active_pos, key=lambda x: x[0])

    def simulate_digest(self, enzymes):
        """Returns lengths in physical order around the circle."""
        active_sites = self.get_active_sites_sorted(enzymes)
        if not active_sites:
            return [float(self.length)]
        lengths = []
        for i in range(len(active_sites)):
            p1 = active_sites[i][0]
            p2 = active_sites[(i + 1) % len(active_sites)][0]
            dist = (p2 - p1) if p2 > p1 else (self.length - p1 + p2)
            lengths.append(round(dist, 4))
        return lengths #note: NOT sorted, follows ordering

    def get_buckets(self, enzymes):
        """Returns (start_pos, length) in physical order."""
        active_sites = self.get_active_sites_sorted(enzymes)
        if not active_sites:
            return [(0.0, float(self.length))]
        lengths = self.simulate_digest(enzymes)
        buckets = []
        for i in range(len(active_sites)):
            buckets.append((active_sites[i][0], lengths[i]))
        return buckets
    
    def format_map(self):
        """Visual representation: Enzyme1/Enzyme2 @ Pos -> length -> ..."""
        if not self.sites:
            return f"--- Blank Plasmid (Length: {self.length}) ---"
        sorted_keys = sorted(self.sites.keys())
        output = []
        for i in range(len(sorted_keys)):
            pos = sorted_keys[i]
            next_pos = sorted_keys[(i + 1) % len(sorted_keys)]
            enzs = "/".join(sorted(self.sites[pos]))
            dist = (next_pos - pos) if next_pos > pos else (self.length - pos + next_pos)
            output.append(f"[{enzs} @ {pos}] ==({round(dist, 1)})==")
        return " ".join(output)

    def get_sig(self):
        """
        Returns a normalized representation of the plasmid.
        """
        # 1. Get all sites sorted by position
        sorted_sites = sorted(self.sites.items())
        if not sorted_sites: return None
        # 2. Extract relative distances and enzyme sets
        # (Distances between sites and the enzymes at those sites)
        n = len(sorted_sites)
        rotations = []
        for i in range(n):
            # Rotate the list so site 'i' is the start
            rotated = sorted_sites[i:] + sorted_sites[:i]
            # Create a representation: ((enzymes), dist_to_next, (enzymes), dist_to_next...)
            # We normalize positions to relative distances
            forward = []
            for j in range(n):
                pos, enz = rotated[j]
                next_pos = rotated[(j + 1) % n][0]
                dist = (next_pos - pos) % self.length
                forward.append((tuple(sorted(enz)), round(dist, 4)))
            # Forward direction
            rotations.append(tuple(forward))
            # Mirror (Reverse) direction
            # We flip the sequence of distances and enzyme sets
            backward = []
            for j in range(n, 0, -1):
                idx = j % n
                prev_idx = (j - 1) % n
                pos, enz = rotated[idx]
                prev_pos = rotated[prev_idx][0]
                dist = (pos - prev_pos) % self.length
                backward.append((tuple(sorted(rotated[idx][1])), round(dist, 4)))
            rotations.append(tuple(backward))
        # The 'signature' is the smallest version
        return min(rotations)

def coverage_possibilities(actual, simulated):
    """
    Finds all ways to partition the 'actual' fragments into the 'simulated' buckets.
    Example: simulated=[5000, 5000], actual=[3000, 2000, 4000, 1000]
    Returns: [[[3000, 2000], [4000, 1000]], [[4000, 1000], [3000, 2000]]]
    """
    from itertools import permutations
    def can_partition(targets, parts):
        if not targets:
            return [[]] if not parts else []
        target = targets[0]
        res = []
        #try all combinations of 'subfragments' that sum to the current 'superfragment'
        for i in range(1, len(parts) + 1):
            for combo in permutations(parts, i):
                if abs(sum(combo) - target) <= 0.1:
                    #if valid, recurse for the remaining targets
                    remaining_parts = list(parts)
                    for x in combo:
                        remaining_parts.remove(x)
                    
                    sub_solutions = can_partition(targets[1:], remaining_parts)
                    for sol in sub_solutions:
                        res.append([list(combo)] + sol)
            #break early if we find a match for this target to save time
            #note: this is for finding a non-ambiguous mapping; remove if multiple valid internal orderings matter)
            #if res: break 
        return res

    #use a set of tuples to distinguish results
    actual_sorted = sorted(actual)
    raw_results = can_partition(simulated, actual)
    unique_results = []
    seen = set()
    for res in raw_results:
        res_tuple = tuple(tuple(group) for group in res)
        if res_tuple not in seen:
            unique_results.append(res)
            seen.add(res_tuple)
    return unique_results

def build_chain(digests):
    #this is weird because it's basically trying to implement the intuition of which maps to try for minimal guesswork
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
            score = (len(new_enzymes) + 1)**2 + len(new_frags) - ((len(common_enzymes) + 3)**1.5)
            overlap = len(set(enzymes) & seen_enzymes)
            if seen_enzymes and overlap == 0: #punishes digests with NO anchoring sites
                #if push is too low, we have to deal with iterating over every base pair (for realistic digests)
                score += 1000  #hard shove to the end, assumes few digests
                #score += (len(digests))**2
            scored.append((score, d))
        scored.sort(key=lambda x: x[0])
        best = scored[0][1]
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
        if not (0.9*average < length < 1.1*average):
            return False
    return average
    '''A brief note on why there's tolerance and this weird check in here, if anyone ever reads this.
        A previous version of this was really bad. At the time, I thought it'd be really easy to make but I neglected the fact I hadn't written any code in the past semester.
        As such, it kinda sucked and I underestimated how long it'd take me to come up with this much more intuitive solution (the unintuitive solution used a bad data structure)
        The second idea was comparing solution spaces for a hopeful O(n) solution with horrible memory complexity.
        Both of these were meant to be integrated with a cnn to use as an android app because I need to learn that
        Also because using the ladder is probably the most annoying thing in b12.
        As such, I predicted that I'd probably have either extensive preprocessing or I could build tolerance into this, before I used question examples from online lectures.
        Some of these functions were just correct so I copied them from my previous file, and they'll come in handy if I ever want to try ml.
    '''

def build_canonical(chain):
    global solutions
    solutions = []
    length = find_length(chain)
    if not length:
        print("Bad digests: Total lengths are inconsistent.")
        return None
    #initialize the plasmid as one continuous 'fragment' with no enzymes yet
    #position starts at 0, length is the total plasmid size
    initial_map = Mapping(length=length)
    result_map = recurse(initial_map, chain)
    unique_maps = []
    seen_signatures = set()
    for mapping in solutions:
        sig = mapping.get_sig() # Ensure this method name matches
        if sig not in seen_signatures:
            unique_maps.append(mapping)
            seen_signatures.add(sig)
    return unique_maps

def recurse(current_mapping, remaining_digests):
    global solutions
    
    # BASE CASE
    if not remaining_digests:
        solutions.append(current_mapping)
        # We return True to indicate a path was found, but the parent
        # loop will keep running to find OTHER paths.
        return True 

    enzymes_in_digest, actual_frags = remaining_digests[0]
    shared_enzymes = set(enzymes_in_digest) & current_mapping.known_enzymes
    new_enzymes = list(set(enzymes_in_digest) - current_mapping.known_enzymes)

    buckets = current_mapping.get_buckets(shared_enzymes)
    simulated_lengths = [b[1] for b in buckets]
    possibilities = coverage_possibilities(actual_frags, simulated_lengths)
    
    if not possibilities:
        return False

    found_at_least_one = False

    # VERIFY/EXPLORE
    for division in possibilities:
        if not new_enzymes:
            if all(len(subgroup) == 1 for subgroup in division):
                if recurse(current_mapping, remaining_digests[1:]):
                    found_at_least_one = True
            continue # Try next division

        num_internal_cuts = sum(len(subgroup) - 1 for subgroup in division)
        if not shared_enzymes:
            num_internal_cuts = len(division[0])

        for enz_setup in product(new_enzymes, repeat=num_internal_cuts):
            trial_map = copy.deepcopy(current_mapping)
            new_sites_to_add = []
            enz_iter = iter(enz_setup)
            
            for i, subgroup in enumerate(division):
                curr_abs_pos, _ = buckets[i]
                pieces = subgroup if not shared_enzymes else subgroup[:-1]
                for length in pieces:
                    curr_abs_pos = (curr_abs_pos + length) % trial_map.length
                    new_sites_to_add.append((curr_abs_pos, next(enz_iter)))
            
            trial_map.add_sites(new_sites_to_add)

            if Counter(trial_map.simulate_digest(enzymes_in_digest)) == Counter(actual_frags):
                # IMPORTANT: We do NOT return immediately. 
                # We record that we found one and let the loop continue.
                if recurse(trial_map, remaining_digests[1:]):
                    found_at_least_one = True
                    if not all_solutions: # If user only wants the first one found
                        return True

    return found_at_least_one

def display_mappings(solution):
    for mapping in solution:
        print(mapping.format_map())

if __name__=="__main__":
    A = (["EcoRI"], [6000, 4000])
    B = (["BamHI"], [5000, 3000, 2000])
    C = (["HindIII"], [7000, 3000])
    D = (["EcoRI", "BamHI"], [3000, 2000, 2000, 1000, 2000])
    E = (["BamHI", "HindIII", "PstI"], [3000, 2000, 2000, 1500, 1500])
    F = (["EcoRI", "HindIII", "XhoI"], [3000, 2500, 2000, 1500, 1000])
    test = [A, B, C, D, E, F]

    #print(coverage_possibilities(E[1], B[1]))

    chain = build_chain(test)
    print("Chain order:", chain)
    print()
    maps = build_canonical(chain)
    display_mappings(maps)
