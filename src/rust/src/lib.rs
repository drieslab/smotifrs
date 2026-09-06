use extendr_api::prelude::*;
use extendr_api::Result;

mod canonical;
mod conditional;
mod enrich;
mod esu;
mod motifclass;
mod graph;
mod motifs;
mod parquet_io;

use crate::graph::SpatialGraphRs;
use crate::motifs::{
    size2::find_size2, size3::find_size3, size4::find_size4, MotifResult,
};

/// Convert `Vec<Option<String>>` to a character vector with `NA`s
/// preserved (R-side reads `NA_character_`).
fn opt_strings_to_robj(v: &[Option<String>]) -> Robj {
    let rstrs: Vec<Rstr> = v.iter().map(|o| match o {
        Some(s) => Rstr::from(s.clone()),
        None    => Rstr::na(),
    }).collect();
    Strings::from_values(rstrs).into()
}

fn motif_result_to_list(r: &MotifResult) -> List {
    let cat = list!(
        motif_id            = r.catalog.motif_id.clone(),
        size                = r.catalog.size.clone(),
        canonical_iso_class = r.catalog.canonical_iso_class.clone(),
        color_tuple         = opt_strings_to_robj(&r.catalog.color_tuple),
        count               = r.catalog.count.clone(),
    );
    let inc = list!(
        motif_id    = r.incidence.motif_id.clone(),
        instance_id = r.incidence.instance_id.clone(),
        cell_id     = r.incidence.cell_id.clone(),
        orbit_id    = r.incidence.orbit_id.clone(),
    );

    let mut imeta_pairs: Vec<(String, Robj)> = Vec::new();
    if let Some(em) = &r.edge_meta {
        if !em.instance_id.is_empty() {
            let v: List = list!(
                instance_id = em.instance_id.clone(),
                v1 = em.v1.clone(),
                v2 = em.v2.clone(),
            );
            imeta_pairs.push(("edge".to_string(), v.into()));
        }
    }
    if let Some(tm) = &r.triangle_meta {
        if !tm.instance_id.is_empty() {
            let v: List = list!(
                instance_id = tm.instance_id.clone(),
                v1 = tm.v1.clone(),
                v2 = tm.v2.clone(),
                v3 = tm.v3.clone(),
            );
            imeta_pairs.push(("triangle".to_string(), v.into()));
        }
    }
    if let Some(wm) = &r.wedge_meta {
        if !wm.instance_id.is_empty() {
            let v: List = list!(
                instance_id = wm.instance_id.clone(),
                end1 = wm.end1.clone(),
                end2 = wm.end2.clone(),
                center = wm.center.clone(),
            );
            imeta_pairs.push(("wedge".to_string(), v.into()));
        }
    }
    let mut size4_classes: Vec<&String> = r.size4_meta.keys().collect();
    size4_classes.sort();
    for cls in size4_classes {
        let m = &r.size4_meta[cls];
        if m.instance_id.is_empty() { continue; }
        let v: List = list!(
            instance_id = m.instance_id.clone(),
            w1 = m.w1.clone(), w2 = m.w2.clone(),
            w3 = m.w3.clone(), w4 = m.w4.clone(),
            d1 = m.d1.clone(), d2 = m.d2.clone(),
            d3 = m.d3.clone(), d4 = m.d4.clone(),
        );
        imeta_pairs.push((cls.clone(), v.into()));
    }
    let imeta = List::from_pairs(imeta_pairs);

    List::from_pairs([
        ("catalog".to_string(),       Robj::from(cat)),
        ("incidence".to_string(),     Robj::from(inc)),
        ("instance_meta".to_string(), Robj::from(imeta)),
    ])
}

fn dispatch_motifs(
    g: &SpatialGraphRs,
    size: i32,
    colored: bool,
    anchored_on: Option<&[u32]>,
) -> std::result::Result<MotifResult, String> {
    match size {
        2 => Ok(find_size2(g, colored, anchored_on)),
        3 => Ok(find_size3(g, colored, anchored_on)),
        4 => Ok(find_size4(g, colored, anchored_on)),
        _ => Err(format!("unsupported size: {}", size)),
    }
}

fn integers_to_u32(ints: &Integers) -> Vec<u32> {
    ints.iter().map(|i| (i.0 - 1) as u32).collect()
}

/// Internal FFI: find motifs from R-supplied vectors. Use the
/// user-facing [`find_motifs_rs()`] wrapper instead.
///
/// @keywords internal
#[extendr]
fn rs_find_motifs(
    cell_ids: Vec<String>,
    cell_types: Vec<String>,
    edge_source: Integers,
    edge_target: Integers,
    size: i32,
    colored: bool,
    anchored_on: Robj,
) -> Result<List> {
    let n = cell_ids.len();
    if cell_types.len() != n {
        return Err(Error::Other(format!(
            "cell_types length {} != cell_ids length {}",
            cell_types.len(), n
        )));
    }
    if edge_source.len() != edge_target.len() {
        return Err(Error::Other("edge_source/edge_target length mismatch".into()));
    }
    let src = integers_to_u32(&edge_source);
    let tgt = integers_to_u32(&edge_target);
    let g = SpatialGraphRs::build(cell_ids, cell_types, &src, &tgt)
        .map_err(Error::Other)?;

    let anchors_vec: Option<Vec<u32>> = if anchored_on.is_null() {
        None
    } else {
        let ints = Integers::try_from(anchored_on)
            .map_err(|e| Error::Other(format!("anchored_on: {:?}", e)))?;
        Some(integers_to_u32(&ints))
    };
    let result = dispatch_motifs(&g, size, colored, anchors_vec.as_deref())
        .map_err(Error::Other)?;
    Ok(motif_result_to_list(&result))
}

/// Internal FFI: read nodes/edges parquet directly in Rust and find
/// motifs. Use the user-facing [`find_motifs_from_parquet()`] wrapper.
///
/// @keywords internal
#[extendr]
fn rs_find_motifs_from_parquet(
    nodes_path: String,
    edges_path: String,
    size: i32,
    colored: bool,
    anchored_on_idx: Robj,
) -> Result<List> {
    let g = parquet_io::graph_from_parquet(&nodes_path, &edges_path)
        .map_err(Error::Other)?;
    let anchors_vec: Option<Vec<u32>> = if anchored_on_idx.is_null() {
        None
    } else {
        let ints = Integers::try_from(anchored_on_idx)
            .map_err(|e| Error::Other(format!("anchored_on_idx: {:?}", e)))?;
        Some(integers_to_u32(&ints))
    };
    let result = dispatch_motifs(&g, size, colored, anchors_vec.as_deref())
        .map_err(Error::Other)?;
    Ok(motif_result_to_list(&result))
}

/// Internal FFI: enumerate motifs and test enrichment in one call, returning
/// only per-class statistics. Per-instance data never crosses into R.
///
/// Cell types arrive as 1-based integer codes plus a level table so that R owns
/// the level ordering; this avoids allocating one string per cell and keeps the
/// two backends from disagreeing under different locale collations.
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_motif_enrichment(
    n_nodes: i32,
    type_codes: Integers,
    type_levels: Vec<String>,
    edge_source: Integers,
    edge_target: Integers,
    size: i32,
    n_perm: i32,
    seed: i32,
    strata: Robj,
    anchored_on: Robj,
    null_kind: &str,
    cond_temp: f64,
) -> Result<List> {
    let n = n_nodes as usize;
    if type_codes.len() != n {
        return Err(Error::Other(format!(
            "type_codes length {} != n_nodes {}",
            type_codes.len(),
            n
        )));
    }
    if edge_source.len() != edge_target.len() {
        return Err(Error::Other("edge_source/edge_target length mismatch".into()));
    }
    if !(2..=4).contains(&size) {
        return Err(Error::Other(format!("size must be 2, 3 or 4; got {}", size)));
    }
    let colors: Vec<u32> = type_codes.iter().map(|i| (i.0 - 1) as u32).collect();
    let src = integers_to_u32(&edge_source);
    let tgt = integers_to_u32(&edge_target);
    let g = SpatialGraphRs::build_coded(n, colors, type_levels, &src, &tgt)
        .map_err(Error::Other)?;

    let strata_vec: Option<Vec<u32>> = if strata.is_null() {
        None
    } else {
        let ints = Integers::try_from(strata)
            .map_err(|e| Error::Other(format!("strata: {:?}", e)))?;
        if ints.len() != n {
            return Err(Error::Other(format!(
                "strata length {} != n_nodes {}",
                ints.len(),
                n
            )));
        }
        Some(ints.iter().map(|i| i.0 as u32).collect())
    };
    let anchors: Option<Vec<u32>> = if anchored_on.is_null() {
        None
    } else {
        let ints = Integers::try_from(anchored_on)
            .map_err(|e| Error::Other(format!("anchored_on: {:?}", e)))?;
        Some(integers_to_u32(&ints))
    };

    let nk = match null_kind {
        "label" => crate::enrich::NullKind::Label,
        "stratified" => crate::enrich::NullKind::Stratified,
        "conditional" => crate::enrich::NullKind::Conditional,
        other => {
            return Err(Error::Other(format!(
                "unknown null '{}'; expected label, stratified or conditional",
                other
            )))
        }
    };
    if nk == crate::enrich::NullKind::Stratified && strata_vec.is_none() {
        return Err(Error::Other(
            "the stratified null needs a strata vector".into(),
        ));
    }

    let r = crate::enrich::run_enrichment(
        &g,
        size as usize,
        n_perm as usize,
        seed as u64,
        strata_vec.as_deref(),
        anchors.as_deref(),
        nk,
        cond_temp,
    )
    .map_err(Error::Other)?;

    let k = size as usize;
    let nclass = r.observed.len();
    // colors flattened row-major, k per class; R reshapes to a matrix
    let mut colors_flat: Vec<Rstr> = Vec::with_capacity(nclass * k);
    let mut motif_id: Vec<String> = Vec::with_capacity(nclass);
    for (i, cc) in r.color_codes.iter().enumerate() {
        let names: Vec<&str> = cc
            .iter()
            .map(|&c| g.color_levels[c as usize].as_str())
            .collect();
        for nm in &names {
            colors_flat.push(Rstr::from(*nm));
        }
        motif_id.push(format!(
            "size{}_{}_{}",
            k,
            r.topo_name[i],
            names.join("-")
        ));
    }

    Ok(list!(
        motif_id = motif_id,
        topology = r.topo_name.clone(),
        size = vec![size; nclass],
        color_flat = Strings::from_values(colors_flat),
        observed = r.observed.clone(),
        mean_null = r.mean_null.clone(),
        sd_null = r.sd_null.clone(),
        p_enrich = r.p_enrich.clone(),
        p_deplete = r.p_deplete.clone(),
        n_instances = r.n_instances as f64,
        n_perm = n_perm,
        k = size,
        cond_dev = r.cond_dev,
        n_edges = r.n_edges as f64,
        cond_accept = r.cond_accept,
        cond_moved = r.cond_moved,
        cond_moves_per_draw = r.cond_moves_per_draw,
        n_nodes_out = n_nodes
    ))
}

/// Internal FFI: vertices of instances belonging to named motif classes.
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_motif_instances(
    n_nodes: i32,
    type_codes: Integers,
    type_levels: Vec<String>,
    edge_source: Integers,
    edge_target: Integers,
    size: i32,
    motif_ids: Vec<String>,
    max_per_class: i32,
) -> Result<List> {
    let n = n_nodes as usize;
    if !(2..=4).contains(&size) {
        return Err(Error::Other(format!("size must be 2, 3 or 4; got {}", size)));
    }
    let colors: Vec<u32> = type_codes.iter().map(|i| (i.0 - 1) as u32).collect();
    let src = integers_to_u32(&edge_source);
    let tgt = integers_to_u32(&edge_target);
    let g = SpatialGraphRs::build_coded(n, colors, type_levels, &src, &tgt)
        .map_err(Error::Other)?;
    let (verts, which) = crate::enrich::collect_instances(
        &g,
        size as usize,
        &motif_ids,
        if max_per_class <= 0 { usize::MAX } else { max_per_class as usize },
    )
    .map_err(Error::Other)?;
    Ok(list!(
        verts = verts,
        which = which.iter().map(|w| *w as i32 + 1).collect::<Vec<i32>>(),
        k = size
    ))
}

/// Internal smoke-test entry point used during package build.
///
/// @keywords internal
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

extendr_module! {
    mod smotifrs;
    fn hello_world;
    fn rs_find_motifs;
    fn rs_find_motifs_from_parquet;
    fn rs_motif_enrichment;
    fn rs_motif_instances;
}
