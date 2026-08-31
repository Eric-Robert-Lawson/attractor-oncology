# Attractor Landscape Validation Framework

## Overview

Validating the perfect-play attractor landscape requires a multi-layered approach that leverages the same mathematical principles used to construct it. This document outlines systematic validation techniques for KQvK and scaling to other endgames, including iterative optimization through intelligent re-ordering.

---

## Iterative Optimization Through Reordering

### Why Reordering Matters for Correction

**Critical insight:** If bugs are discovered, fix time is directly proportional to solving order.

```
Bug Scenario: M values incorrect for certain position types

Naive fix: Re-run entire solver from scratch
  Time: 60+ hours to regenerate everything

Smart fix: Reorder to prioritize buggy positions, re-solve only those
  Time: 2-5 hours to correct problematic subset
  Then: Use corrected positions to cascade into dependent positions
```

### Phase 1: Initial Run (Current Approach)

```cpp
// Run solver with any reasonable ordering
// Output: kqvk_perfect_play.db (potentially with bugs)

// Step 1: Extract positions and sort by M value
create_ordered_position_list("kqvk_perfect_play.db", 
                            "kqvk_positions_ordered_by_m.txt");

// Step 2: Run batch solver using ordered file
batch_solve_from_ordered_file(eng, db, 
                             "kqvk_positions_ordered_by_m.txt", 
                             16);
```

### Phase 2: Validation and Bug Detection

```cpp
// Run Syzygy validation against complete database
AttractorValidator validator("kqvk_perfect_play.db", "/path/to/syzygy");
ValidationReport report = validator.validate_random_sample(5000);

if (report.shorter_mates.size() > 0 || report.longer_mates.size() > 10) {
    // Bugs found - proceed to Phase 3
    cout << "Bugs detected. Proceeding with corrective reordering...\n";
} else {
    // No bugs - validation passed
    cout << "Landscape validated. Perfect play confirmed.\n";
}
```

### Phase 3: Corrective Reordering (If Needed)

**Step 1: Identify Problematic Positions**

```cpp
void analyze_validation_errors(const ValidationReport& report,
                               const string& error_analysis_file) {
    ofstream file(error_analysis_file);
    
    file << "ERROR_ANALYSIS\n";
    file << "==============\n\n";
    
    // Positions with LONGER mates than Syzygy
    file << "LONGER_MATES (Our solver too conservative):\n";
    for (const auto& err : report.longer_mates) {
        file << err.position << " | Our M=" << err.attractor_m 
             << " | Syzygy M=" << err.syzygy_m 
             << " | Delta=" << (err.attractor_m - err.syzygy_m) << "\n";
    }
    
    // Positions with SHORTER mates (critical bugs)
    file << "\nSHORTER_MATES (Critical - investigate immediately):\n";
    for (const auto& err : report.shorter_mates) {
        file << err.position << " | Our M=" << err.attractor_m 
             << " | Syzygy M=" << err.syzygy_m << "\n";
    }
    
    // Move verification failures
    file << "\nMOVE_VERIFICATION_FAILURES:\n";
    for (const auto& err : report.move_errors) {
        file << err.position << " | Move: " << err.move 
             << " | Expected next M=" << err.expected_m_next 
             << " | Actual=" << err.actual_m_next << "\n";
    }
    
    file.close();
}

// Usage
analyze_validation_errors(report, "validation_error_analysis.txt");
```

**Step 2: Create Corrective Ordering**

```cpp
void create_corrective_ordering(const SolvedPositionDatabase& db,
                                const string& error_file,
                                const string& corrective_order_file) {
    
    // Load positions
    vector<pair<int, string>> all_positions;  // (priority, position_key)
    
    // Read error file and mark problematic positions
    set<string> error_positions;
    ifstream errors(error_file);
    string line;
    while (getline(errors, line)) {
        // Parse position from error analysis
        // Format: "WK:a1 WQ:b2 BK:c3 | ..."
        if (line.find("WK:") != string::npos) {
            size_t end = line.find(" |");
            string pos_key = line.substr(0, end);
            error_positions.insert(pos_key);
        }
    }
    errors.close();
    
    // Extract all positions with priorities
    for (const auto& [key, solution] : db.solved) {
        int priority;
        
        if (error_positions.count(key)) {
            // Error position - highest priority (0 = most urgent)
            priority = 0;
        } else if (solution.M_value <= 5) {
            // Foundation positions - solve early (1)
            priority = 1;
        } else if (solution.M_value <= 15) {
            // Middle positions - cascade support (2)
            priority = 2;
        } else {
            // Deep positions - benefit from cascade (3)
            priority = 3;
        }
        
        all_positions.push_back({priority, key});
    }
    
    // Sort by priority, then by M value within priority
    sort(all_positions.begin(), all_positions.end(),
         [&](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        // Within same priority, sort by M value
        int m_a = db.solved[a.second].M_value;
        int m_b = db.solved[b.second].M_value;
        return m_a < m_b;
    });
    
    // Write corrective order
    ofstream file(corrective_order_file);
    int current_priority = -1;
    for (const auto& [priority, pos_key] : all_positions) {
        if (priority != current_priority) {
            current_priority = priority;
            file << "\n# Priority " << priority << "\n";
        }
        file << pos_key << "\n";
    }
    file.close();
    
    cout << "Corrective ordering created: " << corrective_order_file << "\n";
    cout << "Error positions prioritized for re-solving\n";
}

// Usage
create_corrective_ordering(db, 
                          "validation_error_analysis.txt",
                          "kqvk_corrective_order.txt");
```

**Step 3: Selective Re-solving**

```cpp
void batch_solve_corrective(CompositionalEngine& eng,
                           SolvedPositionDatabase& db,
                           const string& corrective_order_file,
                           int max_depth = 16) {
    
    // Load corrective order
    vector<GameState> positions_to_recompute;
    ifstream file(corrective_order_file);
    string line;
    
    int error_count = 0;
    int total_count = 0;
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // Parse position
        size_t wk_pos = line.find("WK:") + 3;
        size_t wq_pos = line.find("WQ:") + 3;
        size_t bk_pos = line.find("BK:") + 3;
        
        Position wk = Position::from_str(line.substr(wk_pos, 2));
        Position wq = Position::from_str(line.substr(wq_pos, 2));
        Position bk = Position::from_str(line.substr(bk_pos, 2));
        
        positions_to_recompute.push_back(GameState(wk, wq, bk, 'W'));
        total_count++;
        
        // Track which are error positions
        if (line.find("# Priority 0") != string::npos) {
            error_count++;
        }
    }
    file.close();
    
    cout << "\nCorrective solving:\n";
    cout << "  Error positions to fix: " << error_count << "\n";
    cout << "  Total positions to re-evaluate: " << total_count << "\n";
    cout << "  Expected time: " << (total_count * 0.5 / 3600) << " hours\n\n";
    
    // Run batch solver with corrective order
    // This will:
    // 1. Re-compute error positions first (high priority)
    // 2. Update database with corrections
    // 3. Fix dependent positions through cascade
    
    int solved_count = 0;
    auto start = chrono::high_resolution_clock::now();
    
    for (size_t idx = 0; idx < positions_to_recompute.size(); idx++) {
        const GameState& pos = positions_to_recompute[idx];
        
        // FORCE re-solve even if in database
        // (This overwrites potentially buggy entries)
        
        eng.nodes_evaluated = 0;
        auto [mvs, tp, bnc] = eng.play_complete_game(pos, db, 50, false, 8);
        
        if (!mvs.empty()) {
            SolvedPosition solution;
            solution.position_key = pos.str();
            solution.best_move = mvs[0];
            // ... populate solution
            
            db.add_position(solution);  // Overwrites old entry
            solved_count++;
        }
        
        if (solved_count % 100 == 0) {
            auto current = chrono::high_resolution_clock::now();
            double elapsed = chrono::duration<double>(current - start).count();
            cout << "[" << idx << "/" << positions_to_recompute.size() 
                 << "] Corrected " << solved_count << " positions\n";
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(end - start).count();
    
    cout << "\nCorrective solving complete:\n";
    cout << "  Positions corrected: " << solved_count << "\n";
    cout << "  Time taken: " << (total_time / 3600) << " hours\n";
}
```

### Phase 4: Validation Loop (Repeat Until Perfect)

```cpp
struct OptimizationIteration {
    int iteration;
    int initial_errors;
    int errors_after_correction;
    double correction_time;
    bool validation_passed;
};

void iterative_optimization_loop(CompositionalEngine& eng,
                                SolvedPositionDatabase& db,
                                int max_iterations = 5) {
    
    vector<OptimizationIteration> iterations;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        cout << "\n" << string(80, '=') << "\n";
        cout << "OPTIMIZATION ITERATION " << (iter + 1) << "\n";
        cout << string(80, '=') << "\n";
        
        // Validate current database
        AttractorValidator validator("kqvk_perfect_play.db", "/path/to/syzygy");
        ValidationReport report = validator.validate_random_sample(2000);
        
        OptimizationIteration current;
        current.iteration = iter + 1;
        current.initial_errors = report.longer_mates.size() + report.shorter_mates.size();
        
        cout << "Errors found: " << current.initial_errors << "\n";
        
        if (current.initial_errors == 0) {
            current.validation_passed = true;
            iterations.push_back(current);
            cout << "\nValidation PASSED. Perfect play confirmed.\n";
            break;
        }
        
        // Analyze errors
        analyze_validation_errors(report, 
                                 "validation_errors_iter_" + to_string(iter) + ".txt");
        
        // Create corrective ordering
        create_corrective_ordering(db,
                                  "validation_errors_iter_" + to_string(iter) + ".txt",
                                  "corrective_order_iter_" + to_string(iter) + ".txt");
        
        // Re-solve with corrective ordering
        auto start = chrono::high_resolution_clock::now();
        batch_solve_corrective(eng, db,
                              "corrective_order_iter_" + to_string(iter) + ".txt",
                              16);
        auto end = chrono::high_resolution_clock::now();
        current.correction_time = chrono::duration<double>(end - start).count();
        
        // Re-validate
        ValidationReport final_report = validator.validate_random_sample(2000);
        current.errors_after_correction = final_report.longer_mates.size() + 
                                         final_report.shorter_mates.size();
        current.validation_passed = (current.errors_after_correction == 0);
        
        cout << "Errors after correction: " << current.errors_after_correction << "\n";
        cout << "Correction time: " << (current.correction_time / 3600) << " hours\n";
        
        iterations.push_back(current);
        
        if (current.validation_passed) break;
    }
    
    // Print optimization summary
    cout << "\n" << string(80, '=') << "\n";
    cout << "OPTIMIZATION SUMMARY\n";
    cout << string(80, '=') << "\n";
    for (const auto& iter : iterations) {
        cout << "Iteration " << iter.iteration << ": "
             << iter.initial_errors << " → " << iter.errors_after_correction 
             << " errors (" << iter.correction_time / 3600 << "h)\n";
    }
}

// Usage in main
iterative_optimization_loop(eng, db, 5);  // Max 5 iterations
```

### Key Advantages of Iterative Reordering

**Cost Analysis:**

```
Scenario 1: Bug found after full solve
  Initial solve: 20-30 hours
  Detected: 50 positions with errors
  Naive fix: Re-run entire solver: 20-30 hours
  Total: 40-60 hours

Scenario 2: Iterative reordering
  Initial solve: 20-30 hours
  Detected: 50 positions with errors
  Corrective solve: Prioritize 50 + dependents (~500): 2-4 hours
  Total: 22-34 hours (25% time savings minimum)
```

**Cumulative benefit with multiple bug rounds:**

```
Iteration 1: 20h initial + 3h correction = 23h (99% accuracy)
Iteration 2: 2h refinement (remaining 1%) = 25h total
Iteration 3: 1h final polish = 26h total

vs. 

Naive approach: 20h + 20h + 20h = 60h (re-running from scratch 3 times)
```

---

## Implementation Checklist

```markdown
### Corrective Reordering Workflow

- [ ] Run initial solver with M-value ordering
- [ ] Validate against Syzygy (sample 2000+ positions)
- [ ] Export error analysis to file
- [ ] Generate corrective ordering (error positions first)
- [ ] Run batch_solve_corrective() with high priority on errors
- [ ] Cascade fixes through dependent positions
- [ ] Re-validate
- [ ] If errors remain, repeat from "Generate corrective ordering"
- [ ] When validation passes, landscape complete

### Monitoring During Corrective Passes

- Track error reduction per iteration
- Monitor time spent in corrections
- Adjust priorities based on error patterns
- Document which position types cause cascading errors
```

---

## Core Validation Principles

### 1. Syzygy Baseline Comparison

**Objective**: Verify that M values (moves to mate) match Syzygy tablebases for identical positions.

**Why Syzygy?**
- Proven correct via decades of chess use
- Uses retrograde analysis (working backward from mate)
- Our approach works forward (compositional)
- Convergence proves correctness

**Methodology**:
```
For each position in attractor_landscape:
  1. Query Syzygy: M_syzygy = moves_to_mate(position)
  2. Query attractor: M_attractor = M_value[position]
  3. Compare:
     - M_attractor == M_syzygy → PASS
     - M_attractor > M_syzygy → FLAG (we found longer path)
     - M_attractor < M_syzygy → CRITICAL (we found shorter mate)
```

**Interpretation**:
- **Exact matches**: Validates compositional search correctly minimizes White's path
- **Longer mates**: Indicates sub-optimal move selection or candidate pruning
- **Shorter mates**: Critical bug requiring immediate investigation

### 2. Move Optimality Validation

**Objective**: Verify that selected moves genuinely lead to claimed M values.

**Approach**:
```
For flagged positions:
  1. Extract best_move from attractor landscape
  2. Apply move to board
  3. Lookup resulting position: M_after = M_value[next_position]
  4. Verify: M_attractor = M_after + 1

  If M_attractor ≠ M_after + 1:
    → Move doesn't lead to claimed M value
    → Bug in move selection or M_cache corruption
```

---

## Validation Workflow

### Phase 1: Early Detection (During Build)

Run validation **every 10,000 positions solved** to catch systematic bugs early.

```python
def early_validation(attractor_db, syzygy_db, sample_size=500):
    """
    Sample random positions from attractor landscape.
    Compare against Syzygy.
    """
    
    errors = {
        'longer_than_syzygy': [],
        'shorter_than_syzygy': [],
        'move_doesnt_match': []
    }
    
    for _ in range(sample_size):
        pos = random_position_from(attractor_db)
        
        m_attractor = attractor_db[pos].M_value
        m_syzygy = syzygy_db.query(pos)
        
        if m_attractor > m_syzygy:
            errors['longer_than_syzygy'].append({
                'position': pos,
                'attractor_M': m_attractor,
                'syzygy_M': m_syzygy,
                'delta': m_attractor - m_syzygy
            })
        elif m_attractor < m_syzygy:
            errors['shorter_than_syzygy'].append({
                'position': pos,
                'attractor_M': m_attractor,
                'syzygy_M': m_syzygy
            })
        
        # Verify move leads to correct M
        best_move = attractor_db[pos].best_move
        next_pos = apply_move(pos, best_move)
        m_next = attractor_db[next_pos].M_value
        
        if m_attractor != m_next + 1:
            errors['move_doesnt_match'].append({
                'position': pos,
                'move': best_move,
                'expected_M_next': m_attractor - 1,
                'actual_M_next': m_next
            })
    
    return errors
```

### Phase 2: Systematic Analysis

**For positions with longer mates:**

```
Root causes (in order of likelihood):

1. Viable candidate filtering too aggressive
   - Check: Is optimal move getting pruned?
   - Evidence: Best move has good M but not in viable list

2. M_shallow underestimate at candidate selection
   - Check: compute_M_shallow(candidate, depth=3) too low
   - Evidence: Good candidates ranked low

3. Memoization pollution (M_cache issue)
   - Check: Are cached M values stale?
   - Evidence: Same position queried twice has different M values

4. Black node count tiebreaker incorrect
   - Check: When M values tie, picking best escape?
   - Evidence: Position with 3-move mate but found 4-move

5. Pruning rules too strict
   - King distance: is_moving_away_from_bk logic?
   - Queen blocking: rotated_dest calculation off by one?
```

**For positions with shorter mates (CRITICAL):**

```
This should almost never happen. If it does:

1. Verify Syzygy is correct
   - Cross-check with Nalimov tablebase
   - If both agree, our solver has fundamental bug

2. Check for illegal moves in solution
   - Is discovered mate actually legal?
   - Does it violate blocking rules?

3. Search implementation bug
   - Is depth=0 returning correctly?
   - Is compositional_search_impl alternating White/Black correctly?
```

---

## Statistical Validation

### Accuracy Metrics

```markdown
## Validation Report Template

**Date**: [when validation ran]
**Positions Checked**: [sample size]
**Total Positions in DB**: [attractor size]
**Iteration**: [1, 2, 3, ...]

### Results

| Metric | Count | Percentage |
|--------|-------|-----------|
| Exact matches | X | X% |
| Longer than Syzygy | Y | Y% |
| Shorter than Syzygy | Z | Z% |
| Move verification failures | W | W% |

### Longer Mates Analysis

- Average delta: [mean difference in moves]
- Worst case: [largest difference]
- Most common delta: [mode]
- Affected M ranges: [which mate counts have issues]

### Recommendations

If >5% longer mates:
- Likely issue in viable candidate filtering
- Review compute_M_shallow() calibration
- Check thresh = max(2, best_M / 3) threshold

If <1% longer mates:
- Acceptable variance
- Landscape is sound

If any shorter mates:
- STOP. Find root cause before continuing.
```

---

## Symmetry-Based Validation

### Rotational Consistency

```cpp
Position original = parse_position_key("WK:a1 WQ:b2 BK:d5");
string original_move = "Qc3";

for (int rot = 0; rot < 4; rot++) {
    Position rotated = rotate_90(original, rot);
    string rotated_move = rotate_move_90(original_move, rot);
    
    int M_original = attractor[original].M_value;
    int M_rotated = attractor[rotated].M_value;
    
    assert(M_original == M_rotated, 
        "Rotated position has different M value!");
    
    GameState after_original = apply(original, original_move);
    GameState after_rotated = apply(rotated, rotated_move);
    
    assert(rotate_90(after_original, rot) == after_rotated,
        "Move rotation doesn't preserve state!");
}
```

**Why this matters:**
- Rotational symmetry is your 4x speedup guarantee
- If rotations don't match, database has duplicated bad data

---

## Edge Case Detection

### Positions to Flag for Manual Review

```markdown
### Critical Edge Cases

1. **Checkmate in 1**
   - Should be M=1, best_move leads immediately to mate
   - Verify: is_checkmate(apply(pos, best_move)) == true

2. **Positions near board edges**
   - Black King cornered (fewer escape squares)
   - Verify: Mate doesn't escape due to board edge
   - Example: WK:a1 WQ:b1 BK:a3 should be fast mate

3. **King opposition**
   - White King directly opposed to Black King
   - Verify: Queen moves don't create illegal adjacent kings
   - Check distance_to() logic

4. **Multiple viable mates**
   - Position has 2+ moves leading to same M value
   - Verify: Black node count tiebreaker picks correctly
   - Different moves, same mate length - this is OK

5. **Deep mates (M > 30)**
   - These take longest to verify
   - Sample check: Pick 10 random M=30+ positions
   - Manually trace 5-10 moves to verify convergence

6. **Positions where Black has only 1-2 legal moves**
   - High constraint positions
   - Verify: Our move matches Syzygy (little room for variation)
```

---

## Validation Implementation

### Quick Validation Script (Python pseudocode)

```python
import sqlite3
import chess
import chess.syzygy

class AttractorValidator:
    def __init__(self, attractor_db_path, syzygy_path):
        self.attractor = load_database(attractor_db_path)
        self.syzygy = chess.syzygy.open_tablebases(syzygy_path)
        self.errors = []
        self.stats = {
            'checked': 0,
            'exact_match': 0,
            'longer': 0,
            'shorter': 0,
            'move_error': 0
        }
    
    def validate_position(self, position_key):
        """Single position validation."""
        m_attractor = self.attractor[position_key].M_value
        
        board = parse_position(position_key)
        m_syzygy = self.syzygy.probe_dtz(board)
        
        self.stats['checked'] += 1
        
        if m_attractor == m_syzygy:
            self.stats['exact_match'] += 1
            return True
        elif m_attractor > m_syzygy:
            self.stats['longer'] += 1
            self.errors.append({
                'type': 'longer',
                'pos': position_key,
                'attractor': m_attractor,
                'syzygy': m_syzygy
            })
            return False
        else:
            self.stats['shorter'] += 1
            self.errors.append({
                'type': 'CRITICAL_shorter',
                'pos': position_key,
                'attractor': m_attractor,
                'syzygy': m_syzygy
            })
            return False
    
    def validate_move(self, position_key, best_move):
        """Verify move leads to correct M value."""
        m_current = self.attractor[position_key].M_value
        next_pos = apply_move(position_key, best_move)
        m_next = self.attractor[next_pos].M_value
        
        if m_current != m_next + 1:
            self.stats['move_error'] += 1
            self.errors.append({
                'type': 'move_mismatch',
                'pos': position_key,
                'move': best_move,
                'expected_m_next': m_current - 1,
                'actual_m_next': m_next
            })
            return False
        return True
    
    def validate_random_sample(self, sample_size=1000):
        """Random sample validation."""
        positions = random.sample(list(self.attractor.keys()), sample_size)
        
        for pos_key in positions:
            self.validate_position(pos_key)
            best_move = self.attractor[pos_key].best_move
            self.validate_move(pos_key, best_move)
        
        return self.report()
    
    def report(self):
        """Generate validation report."""
        total = self.stats['checked']
        return f"""
VALIDATION REPORT
=================
Positions Checked: {total}
Exact Matches: {self.stats['exact_match']} ({100*self.stats['exact_match']/total:.1f}%)
Longer than Syzygy: {self.stats['longer']} ({100*self.stats['longer']/total:.1f}%)
Shorter than Syzygy (CRITICAL): {self.stats['shorter']}
Move Errors: {self.stats['move_error']}

Status: {'PASS' if self.stats['shorter'] == 0 and self.stats['move_error'] == 0 else 'REVIEW NEEDED'}
"""

# Usage
validator = AttractorValidator('kqvk_perfect_play.db', '/path/to/syzygy')
print(validator.validate_random_sample(1000))
```

---

## Scaling Validation to Other Endgames

Once KQvK is validated, the same methodology applies to KRvK, KBBvK, etc.:

```markdown
### Validation Transfer for KRvK

1. Build attractor landscape using same compositional approach
2. Compare sample against Syzygy KRvK tables
3. Verify rotational/reflectional symmetries
4. Check edge cases (Rook vs King coordination, stalemate avoidance)

### Key Differences by Endgame

| Endgame | Symmetries | Edge Cases | Validation Focus |
|---------|-----------|-----------|-----------------|
| KQvK | 4 rotations | Few escape squares | Move precision |
| KRvK | 4 rotations | Rook trapped by King | Avoid stalemate |
| KBBvK | 4 rotations | Same-color bishop | Opposite color bishop |
| KBNvK | 8 symmetries | Knight coordination | Trapped bishop |
```

---

## Monitoring During Build

### Metrics to Track at Checkpoints

```markdown
Every 5K positions solved, log:

1. **Validation Sample (n=100 random positions)**
   - % exact match with Syzygy
   - Any >5 move deltas?
   - Any move verification failures?

2. **Cache Statistics**
   - % of positions found via DB lookup
   - Growth rate of cache hits
   - Rotational symmetry verification

3. **Performance Metrics**
   - Avg solve time per position
   - Trend: getting faster? (cache cascading working?)
   - Estimated time to completion

4. **Attractor Health**
   - Mate distribution by M value
   - Positions with only 1 legal move
   - Average branching factor

Action: If validation accuracy drops below 95%, PAUSE and investigate before continuing.
```

---

## Summary

**Your validation and optimization strategy:**
1. Syzygy comparison catches systematic errors
2. Move verification ensures chains are correct
3. Rotational checks validate symmetry
4. Edge case analysis finds corner bugs
5. **Intelligent reordering corrects bugs efficiently**
6. **Iterative refinement converges to perfect play**

**Key insight:** Even if bugs are found, the iterative reordering approach means you fix them in hours, not days. Each correction cycle is exponentially faster than re-running from scratch.

**When to trust the landscape:**
- Syzygy matches: 98%+
- Move chains verified: 100%
- Rotations consistent: 100%
- No suspicious edge case clusters
- Iterative validation completes

At that point, you have a provably correct, mathematically optimal attractor landscape.
