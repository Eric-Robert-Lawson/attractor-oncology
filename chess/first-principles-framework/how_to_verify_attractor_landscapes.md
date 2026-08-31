
---

## Statistical Validation

### Accuracy Metrics

```markdown
## Validation Report Template

**Date**: [when validation ran]
**Positions Checked**: [sample size]
**Total Positions in DB**: [attractor size]

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
- Acceptable variance (human error in Syzygy lookup)
- Landscape is sound

If any shorter mates:
- STOP. Find root cause before continuing.
```

---

## Symmetry-Based Validation

### Rotational Consistency

Since you generate rotations automatically, verify they're correct:

```cpp
// Pseudo-code
Position original = parse_position_key("WK:a1 WQ:b2 BK:d5");
string original_move = "Qc3";

for (int rot = 0; rot < 4; rot++) {
    Position rotated = rotate_90(original, rot);
    string rotated_move = rotate_move_90(original_move, rot);
    
    // These should have identical M values
    int M_original = attractor[original].M_value;
    int M_rotated = attractor[rotated].M_value;
    
    assert(M_original == M_rotated, 
        "Rotated position has different M value!");
    
    // And the moves should be related by rotation
    GameState after_original = apply(original, original_move);
    GameState after_rotated = apply(rotated, rotated_move);
    
    assert(rotate_90(after_original, rot) == after_rotated,
        "Move rotation doesn't preserve state!");
}
```

**Why this matters:**
- Rotational symmetry is your 4x speedup guarantee
- If rotations don't match, your database has duplicated bad data

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
| KQvK | 4 rotations, 2 reflections | Few escape squares | Move precision |
| KRvK | 4 rotations, 2 reflections | Rook trapped by King | Avoid stalemate |
| KBBvK | 4 rotations, 2 reflections | Same-color bishop issues | Opposite color bishop |
| KBNvK | 8 symmetries | Knight coordination | Trapped bishop scenarios |
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
   - Positions with only 1 legal move (highly constrained)
   - Average branching factor

Action: If validation accuracy drops below 95%, PAUSE and investigate before continuing.
```

---

## Summary

**Your validation strategy is sound:**
1. Syzygy comparison catches systematic errors
2. Move verification ensures chains are correct
3. Rotational checks validate symmetry
4. Edge case analysis finds corner bugs
5. Statistical monitoring detects drift

**Key insight:** The same principles that construct the landscape (topological measure, geometric symmetry, compositional search) are what validate it. You're not external-validating against unrelated logic—you're checking self-consistency.

**When to trust the landscape:**
- Syzygy matches: 98%+
- Move chains verified: 100%
- Rotations consistent: 100%
- No suspicious edge case clusters

At that point, you have a provably correct attractor landscape.
