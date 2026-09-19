# Examples

Two end-to-end design pipelines. Run `setupPath` first (each script also calls it).

| Script | Architecture | Needs |
|---|---|---|
| [`example01_gridLpvYk.m`](example01_gridLpvYk.m) | Grid-based LPV–YK (Sec. II) | LPVTools, Robust Control Toolbox |
| [`example02_partitionedPolytopicLpvYk.m`](example02_partitionedPolytopicLpvYk.m) | Partitioned polytopic LPV–YK (Sec. III) | YALMIP + SDP solver, Robust Control Toolbox |

Both write a `.mat` design package next to themselves. These are gitignored — regenerate rather than commit them.

## What each script demonstrates

**`example01`** — the decoupling claim. Steps 2 and 3 design the nominal and the five local controllers; step 5 solves the only two LMIs in the scheme. Note that nothing in those LMIs refers to the local controllers, the number of subsets, or the switching signal. That is why a sixth controller could be added later without touching any of it.

The script prints the local γ-performances (Table I in the paper) and their maximum, which is the guaranteed global level.

**`example02`** — the partitioning claim, plus the shared-face identity. Vertex controllers are keyed on their coordinates and designed once per unique vertex, so adjacent subsets provably share the controller on their common face; the script asserts $K_{p,3} = K_{p+1,2}$.

It closes with a sweep of the closed-loop H∞ norm against speed, comparing the partitioned design against a gain-scheduled LPV–YK design over the full hull. Both use identical Youla machinery, so the gap between the two curves isolates the benefit of trajectory-aligned partitioning.

## Expected runtime

`example01` is dominated by five `lpvsyn` calls over a 26-point grid; `example02` by the SDP solves. Both are minutes rather than seconds on a typical laptop. The unit tests in `tests/` run in seconds and need only the Control System Toolbox — use those to check the installation before committing to a full design run.
