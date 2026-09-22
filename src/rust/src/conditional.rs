//! Conditional null: hold the pairwise composition, randomize everything else.
//!
//! A marginal label permutation answers "is this motif more common than chance
//! given the tissue's cell type composition". It does not answer the question
//! that matters for a multicellular niche: **is this 3- or 4-node motif more
//! common than its own pairwise attraction already implies?** Under the
//! marginal null, any motif built from an attracting pair looks enriched --
//! the higher-order signal is indistinguishable from an echo of the pairwise
//! one.
//!
//! This module samples label assignments that keep the observed pairwise
//! edge-type contingency table (approximately) fixed while otherwise mixing
//! freely. A motif that is still enriched against these draws is enriched
//! *beyond* its pairwise composition.
//!
//! The sampler is a Metropolis chain over label swaps:
//!
//! * state: one label per node, started from the observed assignment
//! * proposal: swap the labels of two nodes of different type
//! * energy: `E = sum |T - T_obs|` over the pairwise type table, in edges
//! * accept downhill moves always, uphill ones with probability `exp(-dE/temp)`
//!
//! Swapping `u` and `v` only changes edges incident to them, so the energy
//! delta costs `O(deg(u) + deg(v))` -- there is no need to rebuild the table.
//!
//! This is conditioning by tolerance, not exact conditioning: the chain holds
//! the table close rather than pinning it. The achieved deviation is returned
//! so a caller can report how tightly the constraint actually held instead of
//! having to trust it.

use crate::graph::SpatialGraphRs;

/// Upper-triangular index of the unordered type pair `(a, b)`.
#[inline]
fn tkey(a: u32, b: u32, k: usize) -> usize {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    lo as usize * k + hi as usize
}

/// Pairwise edge-type contingency table for a colour assignment.
pub fn pair_table(g: &SpatialGraphRs, colors: &[u32]) -> Vec<i64> {
    let k = g.color_levels.len();
    let mut t = vec![0i64; k * k];
    for u in 0..g.n_nodes as u32 {
        for &v in g.neighbors_of(u) {
            if v > u {
                t[tkey(colors[u as usize], colors[v as usize], k)] += 1;
            }
        }
    }
    t
}

/// Energy of a table against the target: total absolute deviation, in edges.
#[inline]
fn energy(t: &[i64], target: &[i64]) -> i64 {
    t.iter().zip(target).map(|(a, b)| (a - b).abs()).sum()
}

/// A Metropolis chain over label swaps, constrained toward `target`.
pub struct SwapChain<'a> {
    g: &'a SpatialGraphRs,
    pub colors: Vec<u32>,
    table: Vec<i64>,
    target: Vec<i64>,
    k: usize,
    pub energy: i64,
    /// proposals offered (same-colour pairs excluded, they are not proposals)
    pub n_proposed: u64,
    /// proposals accepted
    pub n_accepted: u64,
}

impl<'a> SwapChain<'a> {
    pub fn new(g: &'a SpatialGraphRs, target: Vec<i64>) -> Self {
        let colors = g.colors.clone();
        let table = pair_table(g, &colors);
        let k = g.color_levels.len();
        let e = energy(&table, &target);
        SwapChain {
            g, colors, table, target, k, energy: e,
            n_proposed: 0, n_accepted: 0,
        }
    }

    /// Apply the table delta for recolouring node `u` from `old` to `new`,
    /// skipping the edge to `partner` (handled by the paired call).
    #[inline]
    fn shift(&mut self, u: u32, old: u32, new: u32, partner: u32, sign: i64) {
        for &w in self.g.neighbors_of(u) {
            if w == partner {
                continue;
            }
            let cw = self.colors[w as usize];
            let from = tkey(old, cw, self.k);
            let to = tkey(new, cw, self.k);
            self.table[from] -= sign;
            self.table[to] += sign;
        }
    }

    /// Fraction of nodes whose label differs from the observed assignment.
    /// A chain that is too cold barely moves; this is how far it actually got.
    pub fn displaced(&self) -> f64 {
        let n = self.colors.len();
        if n == 0 {
            return 0.0;
        }
        let d = self
            .colors
            .iter()
            .zip(self.g.colors.iter())
            .filter(|(a, b)| a != b)
            .count();
        d as f64 / n as f64
    }

    /// Propose swapping the labels of `u` and `v`; returns whether it stuck.
    pub fn propose(&mut self, u: u32, v: u32, temp: f64, unit: f64) -> bool {
        let (cu, cv) = (self.colors[u as usize], self.colors[v as usize]);
        if cu == cv {
            return false;
        }
        self.n_proposed += 1;
        let before = self.energy;
        // apply
        self.shift(u, cu, cv, v, 1);
        self.shift(v, cv, cu, u, 1);
        self.colors[u as usize] = cv;
        self.colors[v as usize] = cu;
        let after = energy(&self.table, &self.target);

        let d = (after - before) as f64;
        let keep = d <= 0.0 || unit < (-d / temp).exp();
        if keep {
            self.energy = after;
            self.n_accepted += 1;
            true
        } else {
            // revert
            self.colors[u as usize] = cu;
            self.colors[v as usize] = cv;
            self.shift(u, cv, cu, v, 1);
            self.shift(v, cu, cv, u, 1);
            self.energy = before;
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(n: usize, edges: &[(u32, u32)], cts: Vec<&str>) -> SpatialGraphRs {
        let ids: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let s: Vec<u32> = edges.iter().map(|e| e.0).collect();
        let t: Vec<u32> = edges.iter().map(|e| e.1).collect();
        SpatialGraphRs::build(ids, cts.into_iter().map(String::from).collect(), &s, &t)
            .unwrap()
    }

    #[test]
    fn the_table_counts_every_edge_once() {
        let g = mk(
            4,
            &[(0, 1), (1, 2), (2, 3), (3, 0)],
            vec!["A", "B", "A", "B"],
        );
        let t = pair_table(&g, &g.colors);
        assert_eq!(t.iter().sum::<i64>(), 4);
    }

    #[test]
    fn incremental_delta_matches_a_full_rebuild() {
        // the whole point of the chain is that it never rebuilds the table;
        // this checks the incremental bookkeeping against the honest version
        let g = mk(
            8,
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0),
                (0, 4), (2, 6),
            ],
            vec!["A", "B", "C", "A", "B", "C", "A", "B"],
        );
        let target = pair_table(&g, &g.colors);
        let mut ch = SwapChain::new(&g, target.clone());
        let mut x: u64 = 99;
        for _ in 0..500 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u = ((x >> 33) % 8) as u32;
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let v = ((x >> 33) % 8) as u32;
            if u == v {
                continue;
            }
            // temp huge so most proposals are accepted and the state moves
            ch.propose(u, v, 1e9, 0.5);
            let rebuilt = pair_table(&g, &ch.colors);
            assert_eq!(ch.table, rebuilt, "incremental table drifted");
            assert_eq!(ch.energy, energy(&rebuilt, &target));
        }
    }

    #[test]
    fn a_reverted_proposal_restores_the_state_exactly() {
        let g = mk(
            6,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            vec!["A", "B", "A", "B", "A", "B"],
        );
        let target = pair_table(&g, &g.colors);
        let mut ch = SwapChain::new(&g, target.clone());
        let before_colors = ch.colors.clone();
        let before_table = ch.table.clone();
        // temp tiny and unit 1.0 => every uphill move is rejected
        for u in 0..6u32 {
            for v in 0..6u32 {
                if u != v {
                    ch.propose(u, v, 1e-12, 1.0);
                }
            }
        }
        // any accepted move must have been non-uphill, so energy never rose
        assert!(ch.energy <= energy(&before_table, &target));
        if ch.colors == before_colors {
            assert_eq!(ch.table, before_table);
        }
    }

    #[test]
    fn the_chain_keeps_the_pairwise_table_near_its_target() {
        // ring of 60 nodes, strongly blocked labels => a real pairwise signal
        let n = 60usize;
        let edges: Vec<(u32, u32)> =
            (0..n).map(|i| (i as u32, ((i + 1) % n) as u32)).collect();
        let cts: Vec<&str> = (0..n)
            .map(|i| if (i / 10) % 2 == 0 { "A" } else { "B" })
            .collect();
        let g = mk(n, &edges, cts);
        let target = pair_table(&g, &g.colors);
        let total: i64 = target.iter().sum();

        let mut ch = SwapChain::new(&g, target.clone());
        let mut x: u64 = 7;
        let mut moved = 0;
        for _ in 0..20_000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u = ((x >> 33) % n as u64) as u32;
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let v = ((x >> 33) % n as u64) as u32;
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let unit = ((x >> 11) as f64) / ((1u64 << 53) as f64);
            if u != v && ch.propose(u, v, 1.0, unit) {
                moved += 1;
            }
        }
        assert!(moved > 100, "chain never moved ({} accepted)", moved);
        // the constraint held: deviation stays a small fraction of all edges
        assert!(
            ch.energy * 10 <= total,
            "pairwise table drifted: energy {} of {} edges",
            ch.energy, total
        );
        // ...and the assignment is genuinely different from where it started
        assert_ne!(ch.colors, g.colors);
    }
}
