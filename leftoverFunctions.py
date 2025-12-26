import copy
from collections import Counter
from itertools import product

tolerance = 0

Digest = tuple[list[str], list[float]]

class f:
    def __init__(self, start=None, end=None, length=0, position=0):
        self.startEnzyme = start
        self.endEnzyme = end
        self.length = length
        self.startPosition = position
    
    def __repr__(self):
        return f"f({self.startEnzyme}, {self.endEnzyme}, {self.length}, {self.startPosition})"

def n_sum(fragments, total):
    groups = []
    n = len(fragments)

    def backtrack(x, remain, current_path):
        if abs(remain) <= tolerance:
            groups.append(current_path.copy())
            return

        if remain < -tolerance:
            return

        for i in range(x, n):
            if remain - fragments[i] >= -tolerance:
                current_path.append(i)
                backtrack(i+1, remain-fragments[i], current_path)
                current_path.pop()
    backtrack(0, total, [])
    return groups

def coverage_possibilities(subfragments, superfragments):
    subfragments = sorted(subfragments, reverse=True)
    superfragments = sorted(superfragments, reverse=True)

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
            backtrack(remaining, superfrag_index + 1, curr + [subset])
    backtrack(subfragments, 0, [])
    return results
                
def find_boundaries(digest1, digest2):
    if tolerance == 0:
        return set(digest1) & set(digest2)

    shared = set()
    sorted1 = sorted(digest1)
    sorted2 = sorted(digest2)
    i = 0
    j = 0

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
        subgroups = n_sum(subfrags, superfrag)
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
            score = (len(new_enzymes) + 1)**2 + len(new_frags) - ((len(common_enzymes) + 3)**1.5)
            overlap = len(set(enzymes) & seen_enzymes)
            if seen_enzymes and overlap == 0:
                score += 1000  # hard shove to the end
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

def build_canonical(chain):
    length = find_length(chain)
    if not length:
        print("Bad digests: Total lengths are inconsistent.")
        return None
    # Initialize the plasmid as one continuous fragment with no enzymes yet
    # Position starts at 0, length is the total plasmid size.
    initial_map = [f(start=None, end=None, length=length, position=0)]
    # Start the recursion
    result_map = recurse(initial_map, chain, chain)
    return result_map

def simulate_digest(current_map, enzymes):
    if not enzymes:
        return [sum(f.length for f in current_map)]

    simulated_frags = []
    current_accumulator = 0.0

    # Walk through fragments and "dissolve" boundaries that aren't in the digest
    for frag in current_map:
        current_accumulator += frag.length
        
        # If the boundary at the END of this fragment is one of the enzymes 
        # we are currently simulating, it's a "real" cut.
        if frag.endEnzyme in enzymes:
            simulated_frags.append(round(current_accumulator, 4))
            current_accumulator = 0.0

    # CIRCULAR WRAP-AROUND:
    # If the last fragment didn't end in an active enzyme, the current_accumulator 
    # contains the 'tail'. This tail must be merged with the 'head' (the first fragment).
    if current_accumulator > 0:
        if not simulated_frags:
            # No cuts at all from these enzymes
            return [sum(f.length for f in current_map)]
        
        # Add the tail to the first fragment (which represents the head)
        simulated_frags[0] = round(simulated_frags[0] + current_accumulator, 4)

    return sorted(simulated_frags)

def get_enzyme_assignments(new_enzymes, num_cuts):
    """
    Generated in-line in recurse using itertools.product.
    This generates every permutation of new enzymes for the internal cuts.
    """
    return list(product(new_enzymes, repeat=num_cuts))

def find_mismatches(current_map_fragments, actual_digest_lengths):
    """
    Compares map fragments against actual experimental lengths.
    Returns: (mismatched_fragments, unused_actual_lengths)
    """
    actual_counts = Counter(actual_digest_lengths)
    mismatched_fragments = []
    
    # Check every fragment currently on the map
    for fragment in current_map_fragments:
        # If the fragment's length exists in the actual digest, 
        # it means this piece was NOT cut further.
        found_match = False
        for length in actual_counts:
            if abs(fragment.length - length) <= tolerance and actual_counts[length] > 0:
                actual_counts[length] -= 1
                found_match = True
                break
        # If it wasn't found, this specific fragment object is a "mismatch" that needs to be subdivided.
        if not found_match:
            mismatched_fragments.append(fragment)

    # The remaining lengths in actual_counts are the "sub-fragments" created by the new enzyme cuts.
    unused_actual_lengths = []
    for length, count in actual_counts.items():
        unused_actual_lengths.extend([length] * count)
    # Validation: The sum of mismatched map fragments must equal the sum of the unused actual lengths (within tolerance).
    if abs(sum(f.length for f in mismatched_fragments) - sum(unused_actual_lengths)) > tolerance:
        return None
    return mismatched_fragments, unused_actual_lengths

def format_map(current_map):
    if not current_map:
        return "Empty Map"
    
    # We build a string showing the start enzyme and the length of each fragment
    parts = []
    for frag in current_map:
        start = frag.startEnzyme if frag.startEnzyme else "???"
        parts.append(f"{start} ={int(frag.length)}=")
    
    # Add the closing enzyme of the last fragment to complete the circle
    final_enz = current_map[-1].endEnzyme if current_map[-1].endEnzyme else "???"
    return " ".join(parts) + f" {final_enz}"


def insert_cuts_multi_enzyme(current_map, division, mismatched_frags, enzyme_setup):
    new_map = []
    enz_iter = iter(enzyme_setup)
    
    # Now replacements correctly maps the IDs of fragments inside the TRIAL map
    replacements = {id(frag): sub_lengths for frag, sub_lengths in zip(mismatched_frags, division)}

    for fragment in current_map:
        # If this fragment isn't one of the targets, just pass it through
        if id(fragment) not in replacements:
            new_map.append(fragment)
            continue
        
        sub_lengths = replacements[id(fragment)]
        current_pos = fragment.startPosition

        # CASE: Seeding the blank plasmid
        if fragment.startEnzyme is None:
            seeding_enz = enzyme_setup[0]
            for length in sub_lengths:
                new_map.append(f(seeding_enz, seeding_enz, length, current_pos))
                current_pos += length
        
        # CASE: Subdividing an existing marked fragment
        else:
            last_new_enz = fragment.startEnzyme
            for i, length in enumerate(sub_lengths):
                start_enz = last_new_enz
                if i == len(sub_lengths) - 1:
                    end_enz = fragment.endEnzyme
                else:
                    try:
                        end_enz = next(enz_iter)
                    except StopIteration:
                        end_enz = enzyme_setup[0]
                    last_new_enz = end_enz
                
                new_map.append(f(start_enz, end_enz, length, current_pos))
                current_pos += length
                
    return new_map

def insert_cuts_multi_group(current_map, division, map_groups_to_cut, enzyme_setup):
    new_map = []
    enz_iter = iter(enzyme_setup)
    
    # map_groups_to_cut: [[frag1, frag2], [frag3], ...]
    # division: [[3000, 2000], [1500, 1500], ...]
    
    # We use the ID of the FIRST fragment in each group as the trigger
    # for the replacement, and mark the others for deletion.
    replacements = {}
    for group, sub_lengths in zip(map_groups_to_cut, division):
        replacements[id(group[0])] = (sub_lengths, group)
        for extra_frag in group[1:]:
            replacements[id(extra_frag)] = (None, None) 

    for fragment in current_map:
        fid = id(fragment)
        if fid not in replacements:
            new_map.append(fragment)
            continue
        
        sub_lengths, group = replacements[fid]
        
        # If this fragment was part of a group but not the leader, skip it (it's "consumed")
        if sub_lengths is None:
            continue
            
        current_pos = fragment.startPosition
        
        # CASE: Initial seeding of the plasmid
        if fragment.startEnzyme is None:
            seeding_enz = enzyme_setup[0]
            for length in sub_lengths:
                new_map.append(f(seeding_enz, seeding_enz, length, current_pos))
                current_pos += length
        else:
            # The 'end enzyme' of our new chain must match the 'end enzyme' 
            # of the LAST fragment in the original group we are replacing.
            final_end_enzyme = group[-1].endEnzyme
            
            last_new_enz = fragment.startEnzyme
            for i, length in enumerate(sub_lengths):
                start_enz = last_new_enz
                if i == len(sub_lengths) - 1:
                    end_enz = final_end_enzyme
                else:
                    try:
                        end_enz = next(enz_iter)
                    except StopIteration:
                        # This happens if a bucket has no new internal cuts
                        end_enz = last_new_enz
                    last_new_enz = end_enz
                
                new_map.append(f(start_enz, end_enz, length, current_pos))
                current_pos += length
                
    return new_map

def find_mismatches_by_length(simulated_lengths, actual_lengths):
    sim_counts = Counter(simulated_lengths)
    act_counts = Counter(actual_lengths)
    
    # Remove fragments that match exactly
    for length in list(sim_counts.keys()):
        if length in act_counts:
            match_count = min(sim_counts[length], act_counts[length])
            sim_counts[length] -= match_count
            act_counts[length] -= match_count

    mismatched_sim = []
    for length, count in sim_counts.items():
        if count > 0: mismatched_sim.extend([length] * count)
        
    mismatched_act = []
    for length, count in act_counts.items():
        if count > 0: mismatched_act.extend([length] * count)

    # SUCCESS: Map matches lab data perfectly
    if not mismatched_sim and not mismatched_act:
        return [], []

    # CONTRADICTION: Lab has pieces, but map has no length to provide them
    if not mismatched_sim and mismatched_act:
        return None

    # CONTRADICTION: Largest lab piece is bigger than largest available map gap
    if mismatched_act and max(mismatched_act) > max(mismatched_sim):
        return None

    # CONTRADICTION: Mass balance failure
    if abs(sum(mismatched_sim) - sum(mismatched_act)) > tolerance:
        return None
        
    return mismatched_sim, mismatched_act

def identify_map_groups_ordered(current_map, simulated_lengths, active_enzymes):
    total_len = sum(f.length for f in current_map)
    
    # Get all start positions of fragments whose START or END is an active enzyme
    cut_positions = []
    for f in current_map:
        if f.startEnzyme in active_enzymes:
            cut_positions.append(round(f.startPosition % total_len, 4))
            
    cut_positions = sorted(list(set(cut_positions)))
    
    if not cut_positions:
        return [current_map]

    physical_buckets = []
    for i, start_pos in enumerate(cut_positions):
        end_pos = cut_positions[(i + 1) % len(cut_positions)]
        target_len = (end_pos - start_pos) if end_pos > start_pos else (total_len - start_pos + end_pos)
        
        # Find the fragment that starts at start_pos
        # Use a small tolerance for float comparison
        found_idx = None
        for idx, f in enumerate(current_map):
            if abs((f.startPosition % total_len) - start_pos) < 0.1:
                found_idx = idx
                break
        
        if found_idx is None:
            return None # This triggers the backtrack safely instead of crashing
            
        group = []
        acc_len = 0
        while acc_len < target_len - 0.1:
            curr_f = current_map[found_idx % len(current_map)]
            group.append(curr_f)
            acc_len += curr_f.length
            found_idx += 1
        physical_buckets.append((round(acc_len, 4), group))

    # Align with simulated_lengths order
    ordered_groups = []
    temp_physical = list(physical_buckets)
    for s_len in simulated_lengths:
        match = next((p for p in temp_physical if abs(p[0] - s_len) < 0.1), None)
        if match:
            ordered_groups.append(match[1])
            temp_physical.remove(match)
        else:
            return None
            
    return ordered_groups

def recurse(current_map, remaining_digests, full_chain):
    input()
    print("\nStarting recurse call")
    print(format_map(current_map))
    print(remaining_digests[0])
    if not remaining_digests:
        return current_map

    enzymes_in_digest, actual_frags = remaining_digests[0]
    known_on_map = {f.startEnzyme for f in current_map if f.startEnzyme} | \
                   {f.endEnzyme for f in current_map if f.endEnzyme}
    
    shared_enzymes = set(enzymes_in_digest) & known_on_map
    new_enzymes = list(set(enzymes_in_digest) - known_on_map)

    # These are our "Super-fragment" lengths
    print(shared_enzymes)
    simulated_lengths = simulate_digest(current_map, shared_enzymes)
    print("Below are simulated lengths for the above enzymes")
    print(simulated_lengths)
    
    # These are all possible ways actual_frags fit into simulated_lengths
    possibilities = coverage_possibilities(actual_frags, simulated_lengths)
    print(possibilities)
    
    # These are the actual fragment objects on the map grouped by length
    map_groups = identify_map_groups_ordered(current_map, simulated_lengths, shared_enzymes)
    print(map_groups)

    if not new_enzymes:
        # If any possibility represents a perfect 1-to-1 match (no sub-cutting)
        if any(all(len(bucket) == 1 for bucket in division) for division in possibilities):
            print(f"Result: [V] Verification successful for {enzymes_in_digest}")
            return recurse(current_map, remaining_digests[1:], full_chain)
        else:
            print(f"Result: [X] Verification failed. No new enzymes, but map doesn't match.")
            return None

    # Only run this if we have map_groups and new_enzymes to place
    if map_groups is None: return None

    for division in possibilities:
        # Calculate cuts: if 1 bucket is split into 3 sub-fragments, we need 2 internal cuts
        num_cuts = sum(len(subgroup) - 1 for subgroup in division)
        if map_groups[0][0].startEnzyme is None:
            num_cuts = len(division[0])

        for enzyme_setup in product(new_enzymes, repeat=num_cuts):
            trial_map = copy.deepcopy(current_map)
            
            # Re-locate objects in the new deepcopy memory
            target_positions = [[f.startPosition for f in g] for g in map_groups]
            frags_in_trial = []
            for pos_list in target_positions:
                frags_in_trial.append([f for f in trial_map if f.startPosition in pos_list])

            # Zip the math (division) with the objects (frags_in_trial)
            updated_map = insert_cuts_multi_group(trial_map, division, frags_in_trial, enzyme_setup)
            
            if updated_map:
                result = recurse(updated_map, remaining_digests[1:], full_chain)
                if result: return result
    return None

# Test data
A = (["EcoRI"], [6000, 4000])
B = (["BamHI"], [5000, 3000, 2000])
C = (["HindIII"], [7000, 3000])
D = (["EcoRI", "BamHI"], [3000, 2000, 2000, 1000, 2000])
E = (["BamHI", "HindIII", "PstI"], [3000, 2000, 2000, 1500, 1500])
F = (["EcoRI", "HindIII", "XhoI"], [3000, 2500, 2000, 1500, 1000])
test = [A, B, C, D, E, F]

print(coverage_possibilities(E[1], B[1]))

chain = build_chain(test)
print("Chain order:", chain)
print()
mapping = build_canonical(chain)
print()
if mapping:
    print("FINAL MAPPING:")
    for frag in mapping:
        print(f"  {frag}")
else:
    print("Something either went wrong OR your digests are invalid.")

#print(coverage_possibilities(D[1], A[1]))
#print(coverage_possibilities(D[1], B[1]))


