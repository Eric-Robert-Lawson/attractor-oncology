# Compositional Topological Search Engine

A novel algorithm for computing perfect play in chess endgames using first-principles topological navigation.

## Compilation

Requires C++17 or later and a modern compiler (g++, clang, etc).

```bash
g++ -std=c++17 -O3 -o compositional_trajectory_solver compositional_trajectory_solver.cpp
```

## Usage

Run without debug output:
```bash
./compositional_trajectory_solver
```

Run with debug output:
```bash
./compositional_trajectory_solver --debug
```

## Performance

- **Python version:** 30 minutes for hard positions
- **C++ version:** ~36 seconds for hard positions (50x speedup)
- Memoization eliminates redundant searches across iterative deepening depths
