import copy
from collections import Counter
from itertools import product
import matplotlib.pyplot as plt
import numpy as np


tol = 0
all_solutions = False
solutions = []

import matplotlib.pyplot as plt
import numpy as np

def draw_plasmid(mapping, chain, current_digest_idx):
    plt.clf()
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_direction(-1)  # Clockwise
    ax.set_theta_offset(np.pi/2) # 0 at top
    
    # Constants for layout
    VIEW_LIMIT = 50      
    INNER_RADIUS = 25    # Pushes rings toward the edge
    RING_WIDTH = 2     # Thinner rings
    GAP = 1           # Spacing between rings
    
    # 1. Draw the "Tracks" for each digest in the chain
    for d_idx, (enzymes, frags) in enumerate(chain):
        r_base = INNER_RADIUS + d_idx * (RING_WIDTH + GAP)
        active_enz = set(enzymes)
        
        if active_enz.issubset(mapping.known_enzymes):
            sim_lengths = mapping.simulate_digest(active_enz)
            active_sites = mapping.get_active_sites_sorted(active_enz)
            
            # --- UPDATED COLORMAP SYNTAX ---
            cmap = plt.get_cmap('tab10')
            color = cmap(d_idx % 10)
            alpha = 1.0 if d_idx == current_digest_idx else 0.2
            
            # Draw the Arcs
            for i, length in enumerate(sim_lengths):
                start_pos = active_sites[i][0]
                start_rad = (start_pos / mapping.length) * 2 * np.pi
                width_rad = (length / mapping.length) * 2 * np.pi
                
                ax.bar(start_rad, RING_WIDTH, width=width_rad, bottom=r_base, 
                       color=color, alpha=alpha, align='edge', edgecolor='none')

            # --- CUT SITE MARKERS ---
            for pos, enz_set in active_sites:
                theta = (pos / mapping.length) * 2 * np.pi
                # Draw a white line across the ring to show the cut
                # We extend it slightly (+0.02) to create a "notch" effect
                ax.vlines(theta, r_base - 0.02, r_base + RING_WIDTH + 0.02, 
                          color='white' if alpha == 1.0 else 'black', 
                          lw=1.2, zorder=10)

    # 2. Draw Site Labels
    outermost_r = INNER_RADIUS + len(chain) * (RING_WIDTH + GAP) + 2
    for pos, enzymes in mapping.sites.items():
        theta = (pos / mapping.length) * 2 * np.pi
        label = "/".join(sorted(enzymes))
        
        # Vertical reference line
        ax.plot([theta, theta], [INNER_RADIUS - 1, outermost_r - 1], 
                color='black', lw=0.5, linestyle=':', alpha=0.3)
        
        ax.text(theta, outermost_r, f"{label}\n{int(pos)}", 
                ha='center', va='center', fontsize=7, 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=0.5))

    # Scaling and Styling
    ax.set_ylim(0, VIEW_LIMIT) 
    ax.set_rticks([])
    ax.set_xticks([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)
    
    plt.title(f"Plasmid Map: {int(mapping.length)} bp\nActive: {chain[current_digest_idx][0]}", pad=25)
    plt.draw()
    plt.pause(0.05)

class Mapping:
    def __init__(self, length = 0):
        self.sites = {}
        self.length = length
        self.known_enzymes = set()
        
    def rotate(self, offset):
        new_sites = {}
        for pos, enz_set in self.sites.items():
            # Move the site and wrap around the length
            new_pos = round((pos + offset) % self.length, 4)
            new_sites[new_pos] = enz_set
        self.sites = new_sites

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

def build_seed(sites, length):
    """
    sites: List of tuples like [ (["HindIII"], 500), (["BamHI", "HindIII"], 4000) ]
    """
    initial = Mapping(length=length)
    for enzymes, position in sites:
        # Utilizing your existing add_site logic which handles sets/strings
        initial.add_site(enzymes, position)
    return initial

def build_canonical(chain, seed_sites=None):
    global solutions
    solutions = []
    
    length = find_length(chain)
    if not length:
        print("Bad digests: Total lengths are inconsistent.")
        return None

    initial_map = Mapping(length=length)
    recurse(initial_map, chain, chain, 0)
    
    unique_maps = []
    seen_signatures = set()
    
    for mapping in solutions:
        final_valid_map = None
        
        if seed_sites:
            # We need to test the map and its mirror image
            # because gel data doesn't tell us 'direction'
            for is_mirrored in [False, True]:
                test_map = copy.deepcopy(mapping)
                if is_mirrored:
                    # Mirror: new_pos = (length - old_pos) % length
                    test_map.sites = {round((length - p) % length, 4): e 
                                      for p, e in test_map.sites.items()}
                
                # Anchor: Try to align the first seed site
                anchor_enzs, anchor_pos = seed_sites[0]
                if isinstance(anchor_enzs, str): anchor_enzs = {anchor_enzs}
                
                # Find all positions of the anchor enzyme in our generated map
                possible_starts = [p for p, e in test_map.sites.items() if anchor_enzs.issubset(e)]
                
                for current_start in possible_starts:
                    # Calculate the rotation shift needed
                    shift = (anchor_pos - current_start) % length
                    
                    # Apply shift to all sites
                    rotated_sites = {round((p + shift) % length, 4): e 
                                     for p, e in test_map.sites.items()}
                    
                    # Check if ALL other seed sites match this rotation
                    all_match = True
                    for s_enzs, s_pos in seed_sites:
                        if isinstance(s_enzs, str): s_enzs = {s_enzs}
                        
                        if not any(abs(p - s_pos) <= 1.5 and s_enzs.issubset(e) 
                                   for p, e in rotated_sites.items()):
                            all_match = False
                            break
                    
                    if all_match:
                        test_map.sites = rotated_sites
                        final_valid_map = test_map
                        break
                if final_valid_map: break
            
            if not final_valid_map:
                continue # No rotation or mirror of this solution matches seeds
            mapping = final_valid_map # Use the aligned version
        
        # Finally, check signature to ensure uniqueness
        sig = mapping.get_sig()
        if sig not in seen_signatures:
            unique_maps.append(mapping)
            seen_signatures.add(sig)
            
    return unique_maps

def recurse(current_mapping, remaining_digests, full_chain, current_depth=0):
    global solutions
    
    # VISUALIZATION:
    # Use min() to prevent IndexError on the final success frame
    viz_idx = min(current_depth, len(full_chain) - 1)
    draw_plasmid(current_mapping, full_chain, viz_idx)
    
    # BASE CASE: All digests placed
    if not remaining_digests:
        solutions.append(current_mapping)
        # Final redraw to show the completed map
        draw_plasmid(current_mapping, full_chain, viz_idx)
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
            # Check if current map satisfies this digest without adding new sites
            if all(len(subgroup) == 1 for subgroup in division):
                if recurse(current_mapping, remaining_digests[1:], full_chain, current_depth + 1):
                    found_at_least_one = True
            continue 

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
                if recurse(trial_map, remaining_digests[1:], full_chain, current_depth + 1):
                    found_at_least_one = True
                    if not all_solutions:
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

    G = (["AluI"], [2200, 4800, 8000])
    H = (["PvuI"], [5900, 9100])
    I = (["HindIII"], [3000, 12000])
    #J = (["SacI"], [2000, 4000, 9000])
    K = (["PvuI", "BamHI"], [1100, 5900, 8000])
    #L = (["XmnI"], [4800, 10200])
    #M = (["SacI", "XmnI"], [300, 500, 4000, 8700])
    N = (["AluI", "BamHI"], [1800, 2200, 3000, 8000])
    O = (["HindIII", "AluI"], [1000, 1200, 2000, 4800, 6000])
    P = (["BamHI"], [15000])
    #A = (["EcoRI"], [4000, 2000])
    #B = (["HindIII"], [3500, 2500])
    #C = (["BamHI"], [3000, 2000, 1000])
    #D = (["EcoRI", "HindIII"], [2000, 2000, 1500, 500])
    #E = (["EcoRI", "BamHI"], [3000, 1000, 1000, 1000])
    #F = (["HindIII", "BamHI"], [2500, 2000, 1000, 500])

    #test = [A, B, C, D, E, F]
    #test = [G, H, I, J, K, L, M, N, O, P]
    test = [G, H, I, K, N, O, P]

    seeds = []
    chain = build_chain(test)
    print("Chain order:", chain)
    print()
    maps = build_canonical(chain, seeds)
    display_mappings(maps)
