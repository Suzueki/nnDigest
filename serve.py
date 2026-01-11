import copy
import hashlib
from flask import Flask, render_template, request, jsonify
from collections import Counter
from itertools import product, permutations

app = Flask(__name__)

class Mapping:
    def __init__(self, length=0):
        self.sites = {}
        self.length = length
        self.known_enzymes = set()

    def add_site(self, enzymes, position):
        if isinstance(enzymes, str): enzymes = {enzymes}
        pos = round(float(position) % self.length, 4)
        if pos not in self.sites: self.sites[pos] = set()
        self.sites[pos].update(enzymes)
        self.known_enzymes.update(enzymes)

    def add_sites(self, cuts):
        for position, enzyme in cuts: self.add_site(enzyme, position)

    def get_active_sites_sorted(self, enzymes):
        active_pos = [(pos, enz_set) for pos, enz_set in self.sites.items() if not enz_set.isdisjoint(enzymes)]
        return sorted(active_pos, key=lambda x: x[0])

    def simulate_digest(self, enzymes):
        active_sites = self.get_active_sites_sorted(enzymes)
        if not active_sites: return [float(self.length)]
        lengths = []
        for i in range(len(active_sites)):
            p1 = active_sites[i][0]
            p2 = active_sites[(i + 1) % len(active_sites)][0]
            dist = (p2 - p1) if p2 > p1 else (self.length - p1 + p2)
            lengths.append(round(dist, 4))
        return lengths

    def get_buckets(self, enzymes):
        active_sites = self.get_active_sites_sorted(enzymes)
        if not active_sites: return [(0.0, float(self.length))]
        lengths = self.simulate_digest(enzymes)
        buckets = []
        for i in range(len(active_sites)):
            buckets.append((active_sites[i][0], lengths[i]))
        return buckets

def get_enz_color(enz_string):
    # Generates a vibrant color based on the enzyme name
    return "#" + hashlib.md5(enz_string.encode()).hexdigest()[:6]

def coverage_possibilities(actual, simulated):
    def can_partition(targets, parts):
        if not targets: return [[]] if not parts else []
        target = targets[0]
        res = []
        for i in range(1, len(parts) + 1):
            for combo in permutations(parts, i):
                if abs(sum(combo) - target) <= 1.0:
                    remaining_parts = list(parts)
                    for x in combo: remaining_parts.remove(x)
                    sub_solutions = can_partition(targets[1:], remaining_parts)
                    for sol in sub_solutions: res.append([list(combo)] + sol)
        return res
    return can_partition(simulated, actual)

def build_chain(digests):
    seen_enzymes, chain, remaining = set(), [], digests.copy()
    while remaining:
        scored = []
        for d in remaining:
            enzymes = d[0]
            new_enz = set(enzymes) - seen_enzymes
            overlap = len(set(enzymes) & seen_enzymes)
            score = (len(new_enz) + 1)**2 - (overlap * 2)
            if seen_enzymes and overlap == 0: score += 1000
            scored.append((score, d))
        scored.sort(key=lambda x: x[0])
        best = scored[0][1]
        seen_enzymes.update(best[0]); chain.append(best); remaining.remove(best)
    return chain

def solve_generator(current_mapping, remaining_digests, full_chain, depth=0):
    viz_idx = min(depth, len(full_chain) - 1)
    sites_data = [{"pos": p, "enz": "/".join(sorted(e)), "color": get_enz_color("/".join(sorted(e)))} 
                  for p, e in current_mapping.sites.items()]
    
    is_solution = len(remaining_digests) == 0
    yield {
        "status": "VALID SOLUTION FOUND" if is_solution else f"Depth {depth}: Testing {full_chain[viz_idx][0]}",
        "sites_data": sites_data,
        "total_length": current_mapping.length,
        "is_solution": is_solution,
        "done": False
    }

    if is_solution: return

    enzymes_in_digest, actual_frags = remaining_digests[0]
    shared_enzymes = set(enzymes_in_digest) & current_mapping.known_enzymes
    new_enzymes = list(set(enzymes_in_digest) - current_mapping.known_enzymes)
    buckets = current_mapping.get_buckets(shared_enzymes)
    possibilities = coverage_possibilities(actual_frags, [b[1] for b in buckets])

    for division in possibilities:
        if not new_enzymes:
            if all(len(sub) == 1 for sub in division):
                yield from solve_generator(current_mapping, remaining_digests[1:], full_chain, depth + 1)
            continue
        num_cuts = sum(len(sub) - 1 for sub in division)
        if not shared_enzymes: num_cuts = len(division[0])
        for enz_setup in product(new_enzymes, repeat=num_cuts):
            trial_map = copy.deepcopy(current_mapping)
            new_sites, enz_iter = [], iter(enz_setup)
            for i, subgroup in enumerate(division):
                curr_pos = buckets[i][0]
                pieces = subgroup if not shared_enzymes else subgroup[:-1]
                for length in pieces:
                    curr_pos = (curr_pos + length) % trial_map.length
                    new_sites.append((curr_pos, next(enz_iter)))
            trial_map.add_sites(new_sites)
            if Counter(trial_map.simulate_digest(enzymes_in_digest)) == Counter(actual_frags):
                yield from solve_generator(trial_map, remaining_digests[1:], full_chain, depth + 1)

active_solvers = {}

@app.route('/')
def index(): return render_template('index.html')

@app.route('/start', methods=['POST'])
def start():
    data = request.json
    raw_digests = [(d['enzymes'].split(), [float(f) for f in d['fragments'].split()]) 
                   for d in data['digests'] if d['enzymes'] and d['fragments']]
    chain = build_chain(raw_digests)
    length = sum(raw_digests[0][1])
    active_solvers[request.remote_addr] = solve_generator(Mapping(length=length), chain, chain, 0)
    return jsonify(next(active_solvers[request.remote_addr]))

@app.route('/next')
def next_step():
    try: return jsonify(next(active_solvers[request.remote_addr]))
    except StopIteration: return jsonify({"status": "Search Complete.", "done": True})

if __name__ == '__main__': app.run(debug=True)