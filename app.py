import uuid
import copy
import hashlib
from flask import Flask, render_template, request, jsonify
from collections import Counter
from itertools import product, permutations

app = Flask(__name__)

# --- MIDDLEWARE FOR SUBDIRECTORY HOSTING ---
class PrefixMiddleware(object):
    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):
        path = environ.get('PATH_INFO', '')
        if path.startswith(self.prefix):
            environ['PATH_INFO'] = path[len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
        return self.app(environ, start_response)

app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/nDigest')

# --- CORE LOGIC ---

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
        return [(active_sites[i][0], lengths[i]) for i in range(len(active_sites))]

def get_enz_color(enz_string):
    return "#" + hashlib.md5(enz_string.encode()).hexdigest()[:6]

def coverage_possibilities(actual, simulated):
    def can_partition(targets, parts):
        if not targets: return [[]] if not parts else []
        target = targets[0]
        res = []
        for i in range(1, len(parts) + 1):
            for combo in permutations(parts, i):
                if abs(sum(combo) - target) <= 1.5:
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
            overlap = len(set(enzymes) & seen_enzymes)
            score = (len(set(enzymes) - seen_enzymes) + 1)**2 - (overlap * 2)
            if seen_enzymes and overlap == 0: score += 1000
            scored.append((score, d))
        scored.sort(key=lambda x: x[0])
        best = scored[0][1]
        seen_enzymes.update(best[0]); chain.append(best); remaining.remove(best)
    return chain

# --- ALIGNMENT HELPER ---

def try_align_to_seeds(mapping, seeds):
    """Attempts to rotate or mirror a mapping to match backbone seeds."""
    if not seeds: return mapping
    length = mapping.length
    
    for is_mirrored in [False, True]:
        test_sites = copy.deepcopy(mapping.sites)
        if is_mirrored:
            test_sites = {round((length - p) % length, 4): e for p, e in test_sites.items()}
        
        # Anchor enzyme and position from the first seed
        anchor_enz, anchor_pos = seeds[0]
        # Find everywhere this enzyme exists in our current solution
        possible_starts = [p for p, e in test_sites.items() if anchor_enz in e]
        
        for start_pos in possible_starts:
            shift = (anchor_pos - start_pos) % length
            rotated = {round((p + shift) % length, 4): e for p, e in test_sites.items()}
            
            # Verify if ALL seeds are satisfied by this rotation
            match_count = 0
            for s_enz, s_pos in seeds:
                if any(abs(p - s_pos) <= 2.0 and s_enz in e for p, e in rotated.items()):
                    match_count += 1
            
            if match_count == len(seeds):
                new_map = Mapping(length=length)
                new_map.sites = rotated
                return new_map
    return None

# --- GENERATOR ---

def solve_generator(current_mapping, remaining_digests, full_chain, depth=0, seeds=None):
    is_solution = len(remaining_digests) == 0
    
    # If it's a solution, try to align it to the seeds (backbone)
    display_mapping = current_mapping
    if is_solution and seeds:
        display_mapping = try_align_to_seeds(current_mapping, seeds)
        if not display_mapping:
            return # This branch cannot satisfy the backbone requirements

    sites_data = [{"pos": p, "enz": "/".join(sorted(e)), "color": get_enz_color("/".join(sorted(e)))}
                  for p, e in display_mapping.sites.items()]

    yield {
        "status": "VALID SOLUTION FOUND" if is_solution else f"Depth {depth}: Testing {full_chain[min(depth, len(full_chain)-1)][0]}",
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
                yield from solve_generator(current_mapping, remaining_digests[1:], full_chain, depth + 1, seeds)
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
                yield from solve_generator(trial_map, remaining_digests[1:], full_chain, depth + 1, seeds)

# --- ROUTES ---

active_solvers = {}


@app.route('/')
def index(): 
    return render_template('index.html')

@app.route('/start', methods=['POST'])
def start():
    data = request.json
    # Generate a unique ID for this specific browser session
    sid = str(uuid.uuid4())
    
    raw_digests = [(d['enzymes'].split(), [float(f) for f in d['fragments'].split()])
                   for d in data['digests'] if d['enzymes'] and d['fragments']]
    
    seed_data = data.get('seeds', [])
    seeds = [(s['enzymes'], float(s['pos'])) for s in seed_data if s.get('enzymes') and s.get('pos')]
    
    chain = build_chain(raw_digests)
    length = sum(raw_digests[0][1])
    
    active_solvers[sid] = solve_generator(Mapping(length=length), chain, chain, 0, seeds)
    
    response = jsonify(next(active_solvers[sid]))
    response.set_cookie('session_id', sid) # Store ID in user's browser
    return response

@app.route('/next')
def next_step():
    sid = request.cookies.get('session_id')
    if not sid or sid not in active_solvers:
        return jsonify({"status": "Session expired. Restarting...", "done": True})
    
    try: 
        return jsonify(next(active_solvers[sid]))
    except (StopIteration, KeyError): 
        return jsonify({"status": "Search Complete.", "done": True})

@app.route('/clear', methods=['POST'])
def clear_session():
    sid = request.cookies.get('session_id')
    if sid in active_solvers:
        del active_solvers[sid]
    return jsonify({"status": "Cleared"})

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)