use extendr_api::prelude::*;
use extendr_api::Result;

mod canonical;
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
}
