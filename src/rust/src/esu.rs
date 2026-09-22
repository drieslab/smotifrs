//! Wernicke's ESU: enumerate every connected induced k-subgraph exactly once.
//!
//! The previous size-4 path enumerated `triangles + wedges`, extended each by
//! one neighbor, and deduplicated the results through an `AHashSet<u128>`. That
//! set is the dominant memory term at scale -- it holds one entry per candidate
//! quad, before any of them can be emitted.
//!
//! ESU emits each connected k-subgraph exactly once *by construction*, so there
//! is no dedupe structure at all. It is also naturally parallel over root
//! vertices and streams its output, which is what lets the caller bound memory
//! with a chunk size instead of materializing everything first.
//!
//! The recursion is Wernicke's, specialized to k <= 4 so the working sets are
//! fixed-size stack arrays:
//!
//! ```text
//! for each root v:
//!     extend({v}, {u in N(v) : u > v}, v)
//!
//! extend(sub, ext, v):
//!     if |sub| == k: emit(sub); return
//!     while ext non-empty:
//!         w = pop(ext)
//!         ext' = ext + {u in N(w) : u > v, u not in sub, u not adjacent to sub}
//!         extend(sub + {w}, ext', v)
//! ```
//!
//! The "not adjacent to sub" test is the exclusive-neighborhood condition; it
//! is what makes each subgraph appear once rather than once per generating
//! order.

use crate::graph::SpatialGraphRs;

/// Enumerate connected induced `k`-subgraphs whose smallest-indexed vertex is
/// `root`, calling `emit` once per subgraph. Restricting to a single root makes
/// the whole enumeration trivially parallel: roots partition the output.
pub fn esu_from_root<F: FnMut(&[u32])>(
    g: &SpatialGraphRs,
    root: u32,
    k: usize,
    emit: &mut F,
) {
    let mut sub = [0u32; 4];
    sub[0] = root;
    // extension candidates: neighbours of the root with a larger index
    let mut ext: Vec<u32> = g
        .neighbors_of(root)
        .iter()
        .copied()
        .filter(|&u| u > root)
        .collect();
    if k == 1 {
        emit(&sub[..1]);
        return;
    }
    extend(g, &mut sub, 1, &mut ext, root, k, emit);
}

fn extend<F: FnMut(&[u32])>(
    g: &SpatialGraphRs,
    sub: &mut [u32; 4],
    n_sub: usize,
    ext: &mut Vec<u32>,
    root: u32,
    k: usize,
    emit: &mut F,
) {
    if n_sub == k {
        emit(&sub[..k]);
        return;
    }
    // Iterate over a snapshot: each w is removed from the set passed to the
    // deeper call, which is what stops a subgraph being generated twice.
    while let Some(w) = ext.pop() {
        let mut next: Vec<u32> = ext.clone();
        for &u in g.neighbors_of(w) {
            if u <= root {
                continue;
            }
            if sub[..n_sub].contains(&u) || u == w {
                continue;
            }
            // exclusive neighbourhood: u must not already touch the subgraph
            let mut adj_to_sub = false;
            for &s in &sub[..n_sub] {
                if g.has_edge(s, u) {
                    adj_to_sub = true;
                    break;
                }
            }
            if adj_to_sub || next.contains(&u) {
                continue;
            }
            next.push(u);
        }
        sub[n_sub] = w;
        extend(g, sub, n_sub + 1, &mut next, root, k, emit);
    }
}

/// Adjacency bitmask of the induced subgraph on `verts`, in the bit layout of
/// [`crate::motifclass::pair_bit`].
#[inline]
pub fn induced_mask(g: &SpatialGraphRs, verts: &[u32]) -> u8 {
    let k = verts.len();
    let mut m = 0u8;
    for i in 0..k {
        for j in (i + 1)..k {
            if g.has_edge(verts[i], verts[j]) {
                m |= 1u8 << crate::motifclass::pair_bit(i, j, k);
            }
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::SpatialGraphRs;

    fn mk(n: usize, edges: &[(u32, u32)]) -> SpatialGraphRs {
        let ids: Vec<String> = (0..n).map(|i| format!("c{}", i)).collect();
        let cts: Vec<String> = (0..n).map(|_| "A".to_string()).collect();
        let s: Vec<u32> = edges.iter().map(|e| e.0).collect();
        let t: Vec<u32> = edges.iter().map(|e| e.1).collect();
        SpatialGraphRs::build(ids, cts, &s, &t).unwrap()
    }

    fn count(g: &SpatialGraphRs, k: usize) -> usize {
        let mut n = 0usize;
        for v in 0..g.n_nodes as u32 {
            esu_from_root(g, v, k, &mut |_s| n += 1);
        }
        n
    }

    fn collect(g: &SpatialGraphRs, k: usize) -> Vec<Vec<u32>> {
        let mut out = Vec::new();
        for v in 0..g.n_nodes as u32 {
            esu_from_root(g, v, k, &mut |s| {
                let mut x = s.to_vec();
                x.sort_unstable();
                out.push(x);
            });
        }
        out
    }

    // brute force: every k-subset that is connected, counted once
    fn brute(g: &SpatialGraphRs, k: usize) -> usize {
        let n = g.n_nodes;
        let mut c = 0;
        let idx: Vec<usize> = (0..n).collect();
        fn combos(pool: &[usize], k: usize, start: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
            if cur.len() == k {
                out.push(cur.clone());
                return;
            }
            for i in start..pool.len() {
                cur.push(pool[i]);
                combos(pool, k, i + 1, cur, out);
                cur.pop();
            }
        }
        let mut all = Vec::new();
        combos(&idx, k, 0, &mut Vec::new(), &mut all);
        for s in all {
            // connectivity by BFS within the subset
            let mut seen = vec![false; s.len()];
            seen[0] = true;
            let mut stack = vec![0usize];
            let mut cnt = 1;
            while let Some(a) = stack.pop() {
                for b in 0..s.len() {
                    if !seen[b] && g.has_edge(s[a] as u32, s[b] as u32) {
                        seen[b] = true;
                        cnt += 1;
                        stack.push(b);
                    }
                }
            }
            if cnt == s.len() {
                c += 1;
            }
        }
        c
    }

    #[test]
    fn triangle_has_one_of_each_size() {
        let g = mk(3, &[(0, 1), (1, 2), (0, 2)]);
        assert_eq!(count(&g, 2), 3);
        assert_eq!(count(&g, 3), 1);
    }

    #[test]
    fn path_graph_counts() {
        // P4: 0-1-2-3 has 3 edges, 2 connected triples, 1 connected quad
        let g = mk(4, &[(0, 1), (1, 2), (2, 3)]);
        assert_eq!(count(&g, 2), 3);
        assert_eq!(count(&g, 3), 2);
        assert_eq!(count(&g, 4), 1);
    }

    #[test]
    fn star_k13_counts() {
        let g = mk(4, &[(0, 1), (0, 2), (0, 3)]);
        assert_eq!(count(&g, 2), 3);
        assert_eq!(count(&g, 3), 3);
        assert_eq!(count(&g, 4), 1);
    }

    #[test]
    fn k4_counts() {
        let g = mk(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
        assert_eq!(count(&g, 2), 6);
        assert_eq!(count(&g, 3), 4);
        assert_eq!(count(&g, 4), 1);
    }

    #[test]
    fn every_subgraph_is_emitted_exactly_once() {
        // deterministic pseudo-random graph
        let mut edges = Vec::new();
        let n = 24usize;
        let mut x: u64 = 12345;
        for i in 0..n {
            for j in (i + 1)..n {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if (x >> 33) % 100 < 18 {
                    edges.push((i as u32, j as u32));
                }
            }
        }
        let g = mk(n, &edges);
        for k in 2..=4 {
            let got = collect(&g, k);
            let mut uniq = got.clone();
            uniq.sort();
            uniq.dedup();
            assert_eq!(got.len(), uniq.len(), "k={} emitted a duplicate", k);
            assert_eq!(got.len(), brute(&g, k), "k={} count mismatch", k);
        }
    }

    #[test]
    fn induced_mask_matches_the_graph() {
        let g = mk(4, &[(0, 1), (1, 2), (2, 3), (0, 3)]);
        let m = induced_mask(&g, &[0, 1, 2, 3]);
        assert_eq!(m.count_ones(), 4);
        let (canon, _) = crate::motifclass::canonical_structure(m, 4);
        let topo = crate::motifclass::topologies(4);
        let t = topo.iter().find(|t| t.canon == canon).unwrap();
        assert_eq!(t.name, "cycle");
    }
}
