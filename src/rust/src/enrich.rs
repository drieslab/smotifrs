//! Enumerate-once, recount-many motif enrichment.
//!
//! Under a label-permutation null the graph topology never changes, so the set
//! of motif *instances* is fixed. This module enumerates instances once, keeps
//! them as flat integer arrays, and makes every permutation a relabel-and-count
//! over that fixed set. No strings are constructed anywhere in the loop, and
//! the per-instance data never crosses into R.
//!
//! Two design points carry most of the performance:
//!
//! * **Structure once, colors many.** Canonicalizing a subgraph's structure
//!   costs up to `k!` permutations and happens once per instance at enumeration
//!   time; the result is stored as a topology id plus the vertices in canonical
//!   slot order. Per permutation, only the colors are re-minimized, over `Aut`,
//!   which has 2-8 elements for every topology except K4.
//! * **No null matrix.** Observed counts are known before the permutation pass,
//!   so each draw is folded straight into running `sum`, `sumsq`, `n_ge` and
//!   `n_le`. Memory is O(classes), not O(classes x permutations) -- the latter
//!   is 2 GB at 50k size-4 classes and 10k draws.
//!
//! Classes are addressed densely: a canonical color tuple in base-C, offset by
//! topology. Only classes actually realized are ever touched, tracked through a
//! touched-list so the per-draw fold is proportional to realized classes rather
//! than to the whole address space.

use crate::conditional::{pair_table, SwapChain};
use crate::esu::{esu_from_root, induced_mask};
use crate::graph::SpatialGraphRs;
use crate::motifclass::{
    canonical_color_key, canonical_structure, topologies, Topology,
};
use ahash::AHashMap;
use rayon::prelude::*;

/// Upper bound on the dense class address space (topologies x C^k). Beyond
/// this the scratch array per worker stops being reasonable; the caller is told
/// to merge rare cell types rather than being handed a multi-GB allocation.
const MAX_DENSE: usize = 1 << 24;

/// Motif instances in canonical slot order.
pub struct InstanceStore {
    pub k: usize,
    /// index into the topology table, one per instance
    pub topo: Vec<u8>,
    /// `k` vertex ids per instance, in canonical slot order
    pub verts: Vec<u32>,
}

impl InstanceStore {
    pub fn len(&self) -> usize {
        self.topo.len()
    }
    pub fn is_empty(&self) -> bool {
        self.topo.is_empty()
    }
}

/// Enumerate all connected induced k-subgraphs, canonicalize each one's
/// structure, and store the vertices in canonical slot order. Parallel over
/// ESU roots, which partition the output.
pub fn build_instances(
    g: &SpatialGraphRs,
    k: usize,
    topos: &[Topology],
    anchors: Option<&[u32]>,
) -> InstanceStore {
    let canon_to_idx: AHashMap<u8, u8> = topos
        .iter()
        .enumerate()
        .map(|(i, t)| (t.canon, i as u8))
        .collect();

    let anchor_mask: Option<Vec<bool>> = anchors.map(|a| {
        let mut m = vec![false; g.n_nodes];
        for &v in a {
            if (v as usize) < g.n_nodes {
                m[v as usize] = true;
            }
        }
        m
    });

    let roots: Vec<u32> = (0..g.n_nodes as u32).collect();
    let parts: Vec<(Vec<u8>, Vec<u32>)> = roots
        .par_iter()
        .fold(
            || (Vec::new(), Vec::new()),
            |mut acc, &v| {
                esu_from_root(g, v, k, &mut |sub: &[u32]| {
                    if let Some(m) = &anchor_mask {
                        if !sub.iter().any(|&u| m[u as usize]) {
                            return;
                        }
                    }
                    let mask = induced_mask(g, sub);
                    let (canon, slots) = canonical_structure(mask, k);
                    let ti = match canon_to_idx.get(&canon) {
                        Some(&t) => t,
                        None => return,
                    };
                    acc.0.push(ti);
                    for i in 0..k {
                        acc.1.push(sub[slots[i] as usize]);
                    }
                });
                acc
            },
        )
        .collect();

    let mut topo = Vec::with_capacity(parts.iter().map(|p| p.0.len()).sum());
    let mut verts = Vec::with_capacity(parts.iter().map(|p| p.1.len()).sum());
    for (t, v) in parts {
        topo.extend_from_slice(&t);
        verts.extend_from_slice(&v);
    }
    InstanceStore { k, topo, verts }
}

/// Per-class running statistics over the permutation draws.
#[derive(Clone, Copy, Default)]
struct Acc {
    sum: f64,
    sumsq: f64,
    n_ge: u32,
    n_le: u32,
    n_touched: u32,
}

/// Deterministic RNG so a draw depends only on (seed, draw index), never on
/// how rayon happens to schedule the work.
struct SplitMix64(u64);

impl SplitMix64 {
    #[inline]
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    #[inline]
    fn below(&mut self, n: usize) -> usize {
        (self.next() % (n as u64)) as usize
    }
}

/// Fisher-Yates over the whole vector, or within each stratum when `strata` is
/// given (labels are then only exchanged between cells sharing a stratum).
fn permute(colors: &mut [u32], strata: Option<&Vec<Vec<u32>>>, rng: &mut SplitMix64) {
    match strata {
        None => {
            for i in (1..colors.len()).rev() {
                let j = rng.below(i + 1);
                colors.swap(i, j);
            }
        }
        Some(groups) => {
            for grp in groups {
                for i in (1..grp.len()).rev() {
                    let j = rng.below(i + 1);
                    colors.swap(grp[i] as usize, grp[j] as usize);
                }
            }
        }
    }
}

pub struct EnrichOut {
    pub topo_name: Vec<String>,
    pub color_codes: Vec<Vec<u32>>,
    pub observed: Vec<f64>,
    pub mean_null: Vec<f64>,
    pub sd_null: Vec<f64>,
    pub p_enrich: Vec<f64>,
    pub p_deplete: Vec<f64>,
    pub n_instances: usize,
    /// Mean absolute deviation of the pairwise type table from the observed
    /// one across draws, in edges. Zero for the marginal nulls, which do not
    /// constrain it; for the conditional null it is how tightly the constraint
    /// actually held and should be reported, not assumed.
    pub cond_dev: f64,
    /// Total edges, so `cond_dev` can be read as a fraction.
    pub n_edges: usize,
    /// Conditional null only: fraction of proposals the chain accepted.
    pub cond_accept: f64,
    /// Conditional null only: mean fraction of nodes whose label differs from
    /// the observed assignment at sampling time. A chain held too cold barely
    /// moves, and then every "null" draw is the observed data -- which reads as
    /// "nothing is significant" rather than as a failure. This is the number
    /// that exposes it.
    pub cond_moved: f64,
    /// Conditional null only: accepted swaps between consecutive draws. Draws
    /// are only decorrelated once this is comparable to the node count; a low
    /// value means consecutive draws are near-copies of each other, and of the
    /// observed data.
    pub cond_moves_per_draw: f64,
}

/// Count every instance's canonical class under one color assignment, writing
/// into `counts` and recording first touches in `touched`.
#[inline]
fn tally(
    inst: &InstanceStore,
    colors: &[u32],
    topos: &[Topology],
    stride: usize,
    n_col: usize,
    counts: &mut [u32],
    touched: &mut Vec<u32>,
) {
    let k = inst.k;
    for i in 0..inst.topo.len() {
        let t = inst.topo[i] as usize;
        let base = i * k;
        let mut slot_colors = [0u32; 4];
        for s in 0..k {
            slot_colors[s] = colors[inst.verts[base + s] as usize];
        }
        let key = canonical_color_key(&slot_colors, &topos[t].auts, k);
        // canonical tuple -> base-C address, offset by topology
        let mut addr = 0usize;
        for s in 0..k {
            let c = ((key >> (48 - 16 * s)) & 0xFFFF) as usize;
            addr = addr * n_col + c;
        }
        let idx = t * stride + addr;
        if counts[idx] == 0 {
            touched.push(idx as u32);
        }
        counts[idx] += 1;
    }
}

/// Fold one draw's tallies into the running per-class accumulators, then reset
/// the scratch counts. Shared by the marginal and conditional passes.
#[inline]
fn fold_draw(
    counts: &mut [u32],
    touched: &mut Vec<u32>,
    obs_counts: &[u32],
    acc: &mut AHashMap<u32, Acc>,
) {
    for &idx in touched.iter() {
        let c = counts[idx as usize];
        let o = obs_counts[idx as usize];
        let e = acc.entry(idx).or_default();
        e.sum += c as f64;
        e.sumsq += (c as f64) * (c as f64);
        e.n_touched += 1;
        if c >= o {
            e.n_ge += 1;
        }
        if c <= o {
            e.n_le += 1;
        }
        counts[idx as usize] = 0;
    }
    touched.clear();
}

/// How the null is generated.
#[derive(Clone, Copy, PartialEq)]
pub enum NullKind {
    /// permute labels over all nodes
    Label,
    /// permute labels only within strata
    Stratified,
    /// hold the pairwise edge-type table, randomize the rest
    Conditional,
}

/// Full enrichment run: enumerate once, then fold `n_perm` draws into
/// per-class statistics.
#[allow(clippy::too_many_arguments)]
pub fn run_enrichment(
    g: &SpatialGraphRs,
    k: usize,
    n_perm: usize,
    seed: u64,
    strata: Option<&[u32]>,
    anchors: Option<&[u32]>,
    null_kind: NullKind,
    cond_temp: f64,
) -> Result<EnrichOut, String> {
    let topos = topologies(k);
    let n_col = g.color_levels.len();
    let stride = n_col.checked_pow(k as u32).ok_or("cell type count overflows")?;
    let dense = topos.len().checked_mul(stride).ok_or("class space overflows")?;
    if dense > MAX_DENSE {
        return Err(format!(
            "class address space is {} ({} cell types, size {}), above the {} limit; \
             merge rare cell types or run a smaller motif size",
            dense, n_col, k, MAX_DENSE
        ));
    }

    let inst = build_instances(g, k, &topos, anchors);
    let n_inst = inst.len();

    // observed
    let mut obs_counts = vec![0u32; dense];
    let mut obs_touched: Vec<u32> = Vec::new();
    tally(&inst, &g.colors, &topos, stride, n_col, &mut obs_counts, &mut obs_touched);

    // group node indices by stratum, if stratifying
    let strata_groups: Option<Vec<Vec<u32>>> = strata.map(|s| {
        let mut m: AHashMap<u32, Vec<u32>> = AHashMap::new();
        for (i, &v) in s.iter().enumerate() {
            m.entry(v).or_default().push(i as u32);
        }
        let mut keys: Vec<u32> = m.keys().copied().collect();
        keys.sort_unstable();
        keys.into_iter().map(|kk| m.remove(&kk).unwrap()).collect()
    });

    let reduce_acc = |mut a: AHashMap<u32, Acc>, b: AHashMap<u32, Acc>| {
        for (kk, v) in b {
            let e = a.entry(kk).or_default();
            e.sum += v.sum;
            e.sumsq += v.sumsq;
            e.n_ge += v.n_ge;
            e.n_le += v.n_le;
            e.n_touched += v.n_touched;
        }
        a
    };

    let (merged, cond_dev, cond_accept, cond_moved): (AHashMap<u32, Acc>, f64, f64, f64) =
    if null_kind == NullKind::Conditional {
        // A Markov chain is sequential, so parallelism comes from running one
        // independent chain per worker rather than from splitting the draws.
        let n_chains = rayon::current_num_threads().max(1).min(n_perm);
        let per_chain = n_perm.div_ceil(n_chains);
        let target = pair_table(g, &g.colors);
        let burn = 20 * g.n_nodes;
        let thin = 2 * g.n_nodes;

        let parts: Vec<(AHashMap<u32, Acc>, f64, usize, u64, u64, f64)> = (0..n_chains)
            .into_par_iter()
            .map(|ci| {
                let mut acc = AHashMap::<u32, Acc>::new();
                let mut counts = vec![0u32; dense];
                let mut touched = Vec::<u32>::new();
                let mut rng =
                    SplitMix64(seed ^ (ci as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15));
                let mut chain = SwapChain::new(g, target.clone());
                let nn = g.n_nodes;
                let mut step = |ch: &mut SwapChain, rng: &mut SplitMix64| {
                    let u = rng.below(nn) as u32;
                    let v = rng.below(nn) as u32;
                    let unit = (rng.next() >> 11) as f64 / (1u64 << 53) as f64;
                    if u != v {
                        ch.propose(u, v, cond_temp, unit);
                    }
                };
                for _ in 0..burn {
                    step(&mut chain, &mut rng);
                }
                let mut dev = 0f64;
                let mut moved = 0f64;
                let start = ci * per_chain;
                let end = ((ci + 1) * per_chain).min(n_perm);
                for _ in start..end {
                    for _ in 0..thin {
                        step(&mut chain, &mut rng);
                    }
                    dev += chain.energy as f64;
                    moved += chain.displaced();
                    tally(&inst, &chain.colors, &topos, stride, n_col, &mut counts, &mut touched);
                    fold_draw(&mut counts, &mut touched, &obs_counts, &mut acc);
                }
                (
                    acc, dev, end.saturating_sub(start),
                    chain.n_proposed, chain.n_accepted, moved,
                )
            })
            .collect();

        let total: f64 = parts.iter().map(|p| p.1).sum();
        let ndraw: usize = parts.iter().map(|p| p.2).sum();
        let prop: u64 = parts.iter().map(|p| p.3).sum();
        let acc_n: u64 = parts.iter().map(|p| p.4).sum();
        let movsum: f64 = parts.iter().map(|p| p.5).sum();
        let m = parts
            .into_iter()
            .map(|p| p.0)
            .fold(AHashMap::new(), reduce_acc);
        let d = if ndraw > 0 { total / ndraw as f64 } else { 0.0 };
        let a = if prop > 0 { acc_n as f64 / prop as f64 } else { 0.0 };
        let mv = if ndraw > 0 { movsum / ndraw as f64 } else { 0.0 };
        (m, d, a, mv)
    } else {
        let m = (0..n_perm)
            .into_par_iter()
            .fold(
                || (AHashMap::<u32, Acc>::new(), vec![0u32; dense], Vec::<u32>::new()),
                |(mut acc, mut counts, mut touched), draw| {
                    let mut rng =
                        SplitMix64(seed ^ (draw as u64).wrapping_mul(0x9E3779B97F4A7C15));
                    let mut cols = g.colors.clone();
                    permute(&mut cols, strata_groups.as_ref(), &mut rng);
                    tally(&inst, &cols, &topos, stride, n_col, &mut counts, &mut touched);
                    fold_draw(&mut counts, &mut touched, &obs_counts, &mut acc);
                    (acc, counts, touched)
                },
            )
            .map(|(acc, _, _)| acc)
            .reduce(AHashMap::new, reduce_acc);
        (m, 0.0, f64::NAN, f64::NAN)
    };

    // report every class seen in the observed graph or in any draw
    let mut all: Vec<u32> = obs_touched.clone();
    all.extend(merged.keys().copied());
    all.sort_unstable();
    all.dedup();

    let np = n_perm as f64;
    let mut out = EnrichOut {
        topo_name: Vec::with_capacity(all.len()),
        color_codes: Vec::with_capacity(all.len()),
        observed: Vec::with_capacity(all.len()),
        mean_null: Vec::with_capacity(all.len()),
        sd_null: Vec::with_capacity(all.len()),
        p_enrich: Vec::with_capacity(all.len()),
        p_deplete: Vec::with_capacity(all.len()),
        n_instances: n_inst,
        cond_dev,
        n_edges: g.edge_keys.len(),
        cond_accept,
        cond_moved,
        cond_moves_per_draw: if cond_accept.is_finite() {
            cond_accept * (2 * g.n_nodes) as f64
        } else {
            f64::NAN
        },
    };

    for idx in all {
        let t = (idx as usize) / stride;
        let mut addr = (idx as usize) % stride;
        let mut cols = vec![0u32; k];
        for s in (0..k).rev() {
            cols[s] = (addr % n_col) as u32;
            addr /= n_col;
        }
        let o = obs_counts[idx as usize] as f64;
        let a = merged.get(&idx).copied().unwrap_or_default();

        // draws in which this class never appeared contribute a count of 0
        let n_zero = np - a.n_touched as f64;
        let mean = a.sum / np;
        let var = (a.sumsq / np - mean * mean).max(0.0);
        // 0 >= observed only when observed is 0; 0 <= observed always
        let n_ge = a.n_ge as f64 + if o == 0.0 { n_zero } else { 0.0 };
        let n_le = a.n_le as f64 + n_zero;

        out.topo_name.push(topos[t].name.to_string());
        out.color_codes.push(cols);
        out.observed.push(o);
        out.mean_null.push(mean);
        out.sd_null.push(var.sqrt());
        // unbiased "+1" estimator, both tails
        out.p_enrich.push((1.0 + n_ge) / (1.0 + np));
        out.p_deplete.push((1.0 + n_le) / (1.0 + np));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(n: usize, edges: &[(u32, u32)], cts: Vec<&str>) -> SpatialGraphRs {
        let ids: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let s: Vec<u32> = edges.iter().map(|e| e.0).collect();
        let t: Vec<u32> = edges.iter().map(|e| e.1).collect();
        SpatialGraphRs::build(ids, cts.into_iter().map(String::from).collect(), &s, &t).unwrap()
    }

    #[test]
    fn observed_counts_sum_to_the_instance_count() {
        let g = mk(
            6,
            &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 3)],
            vec!["A", "B", "A", "B", "A", "B"],
        );
        for k in 2..=4 {
            let r = run_enrichment(&g, k, 50, 1, None, None, NullKind::Label, 1.0).unwrap();
            let tot: f64 = r.observed.iter().sum();
            assert_eq!(tot as usize, r.n_instances, "k={}", k);
        }
    }

    #[test]
    fn every_draw_conserves_the_instance_count() {
        // mean_null summed over classes must equal the number of instances:
        // a permutation moves instances between classes, never creates them
        let g = mk(
            8,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0), (0, 4)],
            vec!["A", "B", "C", "A", "B", "C", "A", "B"],
        );
        for k in 2..=4 {
            let r = run_enrichment(&g, k, 200, 7, None, None, NullKind::Label, 1.0).unwrap();
            let tot: f64 = r.mean_null.iter().sum();
            assert!(
                (tot - r.n_instances as f64).abs() < 1e-9,
                "k={} mean_null sums to {} not {}",
                k, tot, r.n_instances
            );
        }
    }

    #[test]
    fn p_values_are_valid_probabilities() {
        let g = mk(
            7,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0)],
            vec!["A", "A", "B", "B", "C", "C", "A"],
        );
        let r = run_enrichment(&g, 3, 99, 3, None, None, NullKind::Label, 1.0).unwrap();
        for i in 0..r.observed.len() {
            assert!(r.p_enrich[i] > 0.0 && r.p_enrich[i] <= 1.0);
            assert!(r.p_deplete[i] > 0.0 && r.p_deplete[i] <= 1.0);
            assert!(r.p_enrich[i] >= 1.0 / 100.0);
            assert!(r.sd_null[i] >= 0.0);
        }
    }

    #[test]
    fn a_single_cell_type_yields_one_class_per_topology() {
        let g = mk(
            5,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
            vec!["A", "A", "A", "A", "A"],
        );
        let r = run_enrichment(&g, 3, 20, 1, None, None, NullKind::Label, 1.0).unwrap();
        // C5 has 5 wedges and no triangles -> exactly one class
        assert_eq!(r.observed.len(), 1);
        assert_eq!(r.topo_name[0], "open");
        assert_eq!(r.observed[0], 5.0);
        // with one color the null cannot move anything
        assert!((r.mean_null[0] - 5.0).abs() < 1e-9);
        assert!(r.sd_null[0] < 1e-9);
    }

    #[test]
    fn results_are_reproducible_and_seed_dependent() {
        let g = mk(
            6,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            vec!["A", "B", "A", "B", "A", "B"],
        );
        let a = run_enrichment(&g, 3, 50, 42, None, None, NullKind::Label, 1.0).unwrap();
        let b = run_enrichment(&g, 3, 50, 42, None, None, NullKind::Label, 1.0).unwrap();
        let c = run_enrichment(&g, 3, 50, 43, None, None, NullKind::Label, 1.0).unwrap();
        assert_eq!(a.mean_null, b.mean_null);
        assert_ne!(a.mean_null, c.mean_null);
    }

    #[test]
    fn stratified_null_never_moves_labels_across_strata() {
        // two disjoint triangles, each a pure single-color stratum: permuting
        // within strata cannot change anything at all
        let g = mk(
            6,
            &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)],
            vec!["A", "A", "A", "B", "B", "B"],
        );
        let strata = vec![0u32, 0, 0, 1, 1, 1];
        let r = run_enrichment(&g, 3, 50, 5, Some(&strata), None, NullKind::Stratified, 1.0).unwrap();
        for i in 0..r.observed.len() {
            assert!((r.mean_null[i] - r.observed[i]).abs() < 1e-9);
            assert!(r.sd_null[i] < 1e-9);
        }
    }
}

#[cfg(test)]
mod scale {
    use super::*;
    use std::time::Instant;

    /// Grid graph with ~6 neighbours per node, the shape of a spatial network.
    fn grid_graph(n: usize, n_types: usize) -> SpatialGraphRs {
        let side = (n as f64).sqrt().ceil() as usize;
        let ids: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let cts: Vec<String> = (0..n).map(|i| format!("T{}", i % n_types)).collect();
        let (mut s, mut t) = (Vec::new(), Vec::new());
        for i in 0..n {
            let (r, c) = (i / side, i % side);
            for (dr, dc) in [(0i64, 1i64), (1, 0), (1, 1), (1, -1)] {
                let (nr, nc) = (r as i64 + dr, c as i64 + dc);
                if nr < 0 || nc < 0 || nc >= side as i64 {
                    continue;
                }
                let j = (nr as usize) * side + (nc as usize);
                if j < n {
                    s.push(i as u32);
                    t.push(j as u32);
                }
            }
        }
        SpatialGraphRs::build(ids, cts, &s, &t).unwrap()
    }

    #[test]
    #[ignore]
    fn scale_report() {
        for &(n, k, np) in &[
            (50_000usize, 3usize, 1000usize),
            (50_000, 4, 1000),
            (170_000, 3, 1000),
            (170_000, 4, 1000),
        ] {
            let g = grid_graph(n, 12);
            let t0 = Instant::now();
            let topos = topologies(k);
            let inst = build_instances(&g, k, &topos, None);
            let t_enum = t0.elapsed().as_secs_f64();
            let n_inst = inst.len();
            drop(inst);

            let t1 = Instant::now();
            let r = run_enrichment(&g, k, np, 1, None, None, NullKind::Label, 1.0).unwrap();
            let t_all = t1.elapsed().as_secs_f64();
            println!(
                "n={:>7} k={} perms={:>5} | edges={:>9} instances={:>12} classes={:>7} \
                 | enum {:>7.2}s  total {:>8.2}s",
                n, k, np, g.edge_keys.len(), n_inst, r.observed.len(), t_enum, t_all
            );
        }
    }
}
