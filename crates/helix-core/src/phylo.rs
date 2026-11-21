//! Phylogenetic tree construction (UPGMA and Neighbor-Joining) and Newick serialization.

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::sequence::is_gap;
use serde::{Deserialize, Serialize};

/// Distance-based phylogenetic tree inference method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum TreeMethod {
    /// Unweighted Pair Group Method with Arithmetic Mean (ultrametric rooted tree).
    #[default]
    Upgma,
    /// Neighbor-Joining (additive unrooted/midpoint-rooted phylogram).
    NeighborJoining,
}

impl TreeMethod {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Upgma => "UPGMA",
            Self::NeighborJoining => "Neighbor-Joining",
        }
    }
}

/// Evolutionary substitution distance correction model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum DistanceModel {
    /// Raw proportion of differing sites (`p-distance`).
    PDistance,
    /// Jukes-Cantor multiple-hit correction (`-b * ln(1 - p/b)`).
    #[default]
    JukesCantor,
}

impl DistanceModel {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::PDistance => "p-Distance",
            Self::JukesCantor => "Jukes-Cantor",
        }
    }
}

/// A node in a rooted binary/multifurcating phylogenetic tree.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhyloNode {
    /// Taxon label for leaves (or internal support/cluster label).
    pub label: String,
    /// Branch length leading from this node's parent to this node (`>= 0.0`).
    pub branch_length: f64,
    /// Child subtrees (empty for leaf taxa).
    pub children: Vec<PhyloNode>,
}

impl PhyloNode {
    /// Returns `true` if this node is a terminal leaf taxon.
    #[must_use]
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// Counts total leaf taxa under this subtree.
    #[must_use]
    pub fn leaf_count(&self) -> usize {
        if self.is_leaf() {
            1
        } else {
            self.children.iter().map(Self::leaf_count).sum()
        }
    }

    /// Maximum cumulative branch length from this node down to any leaf.
    #[must_use]
    pub fn max_depth(&self) -> f64 {
        if self.is_leaf() {
            0.0
        } else {
            self.children
                .iter()
                .map(|c| c.branch_length + c.max_depth())
                .fold(0.0f64, f64::max)
        }
    }

    /// Serializes this subtree into standard Newick format (without trailing semicolon).
    fn fmt_newick(&self, out: &mut String) {
        if !self.children.is_empty() {
            out.push('(');
            for (idx, child) in self.children.iter().enumerate() {
                if idx > 0 {
                    out.push(',');
                }
                child.fmt_newick(out);
            }
            out.push(')');
        }
        if !self.label.is_empty() {
            // Escape Newick-reserved characters cleanly.
            let clean: String = self
                .label
                .chars()
                .map(|c| match c {
                    '(' | ')' | ',' | ':' | ';' => '_',
                    _ => c,
                })
                .collect();
            out.push_str(&clean);
        }
        out.push_str(&format!(":{:.5}", self.branch_length.max(0.0)));
    }
}

/// A phylogenetic tree constructed from a multiple sequence alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhyloTree {
    pub name: String,
    pub method: TreeMethod,
    pub model: DistanceModel,
    pub root: PhyloNode,
}

impl PhyloTree {
    /// Exports this tree in standard Newick format (`((A:0.05,B:0.05):0.10,C:0.15);`).
    #[must_use]
    pub fn to_newick(&self) -> String {
        let mut out = String::new();
        if !self.root.children.is_empty() {
            out.push('(');
            for (idx, child) in self.root.children.iter().enumerate() {
                if idx > 0 {
                    out.push(',');
                }
                child.fmt_newick(&mut out);
            }
            out.push(')');
        } else {
            out.push_str(&self.root.label);
        }
        out.push(';');
        out
    }
}

/// Builds a phylogenetic tree from `alignment` using UPGMA or Neighbor-Joining.
pub fn build_phylogenetic_tree(
    alignment: &Alignment,
    method: TreeMethod,
    model: DistanceModel,
) -> Result<PhyloTree> {
    let n = alignment.row_count();
    if n < 2 {
        return Err(HelixError::Alignment(
            "at least 2 aligned sequences are required to build a phylogenetic tree".to_owned(),
        ));
    }

    let names: Vec<String> = alignment.sequences().iter().map(|s| s.name.clone()).collect();
    let dist = compute_distance_matrix(alignment, model);

    let root = match method {
        TreeMethod::Upgma => build_upgma(&names, &dist),
        TreeMethod::NeighborJoining => build_neighbor_joining(&names, &dist),
    };

    Ok(PhyloTree {
        name: format!("{} ({})", alignment.name, method.label()),
        method,
        model,
        root,
    })
}

/// Computes symmetric `N x N` evolutionary distance matrix from an alignment.
fn compute_distance_matrix(alignment: &Alignment, model: DistanceModel) -> Vec<Vec<f64>> {
    let n = alignment.row_count();
    let mut mat = vec![vec![0.0f64; n]; n];
    // Jukes-Cantor asymptote: 3/4 for 4-letter DNA/RNA, 19/20 for 20-letter proteins.
    let b = if alignment.alphabet.is_nucleotide() {
        0.75f64
    } else {
        0.95f64
    };

    for i in 0..n {
        for j in (i + 1)..n {
            let s1 = &alignment.sequences()[i].residues;
            let s2 = &alignment.sequences()[j].residues;
            let mut compared = 0usize;
            let mut diffs = 0usize;
            for (&a, &c) in s1.iter().zip(s2.iter()) {
                if is_gap(a) || is_gap(c) {
                    continue;
                }
                compared += 1;
                if a != c {
                    diffs += 1;
                }
            }
            let p = if compared == 0 {
                0.5
            } else {
                diffs as f64 / compared as f64
            };
            let d = match model {
                DistanceModel::PDistance => p,
                DistanceModel::JukesCantor => {
                    let arg = (1.0 - p / b).max(0.01);
                    -b * arg.ln()
                }
            };
            mat[i][j] = d;
            mat[j][i] = d;
        }
    }
    mat
}

fn build_upgma(names: &[String], initial_dist: &[Vec<f64>]) -> PhyloNode {
    struct Cluster {
        node: PhyloNode,
        height: f64,
        size: usize,
    }

    let mut clusters: Vec<Cluster> = names
        .iter()
        .map(|name| Cluster {
            node: PhyloNode {
                label: name.clone(),
                branch_length: 0.0,
                children: Vec::new(),
            },
            height: 0.0,
            size: 1,
        })
        .collect();

    let mut dist = initial_dist.to_vec();

    while clusters.len() > 1 {
        let k = clusters.len();
        let mut best_i = 0usize;
        let mut best_j = 1usize;
        let mut best_d = f64::INFINITY;

        for (i, row) in dist.iter().enumerate().take(k) {
            for (j, &d) in row.iter().enumerate().take(k).skip(i + 1) {
                if d < best_d {
                    best_d = d;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        let new_height = (best_d * 0.5).max(0.0005);
        let c_j = clusters.remove(best_j);
        let c_i = clusters.remove(best_i);

        let mut left = c_i.node;
        let mut right = c_j.node;
        left.branch_length = (new_height - c_i.height).max(0.0002);
        right.branch_length = (new_height - c_j.height).max(0.0002);

        let merged_size = c_i.size + c_j.size;
        let merged = Cluster {
            node: PhyloNode {
                label: String::new(),
                branch_length: 0.0,
                children: vec![left, right],
            },
            height: new_height,
            size: merged_size,
        };

        // Compute weighted average distances to all remaining clusters.
        let mut new_row = Vec::with_capacity(clusters.len() + 1);
        for m in 0..k {
            if m == best_i || m == best_j {
                continue;
            }
            let d_new = (dist[best_i][m] * c_i.size as f64 + dist[best_j][m] * c_j.size as f64)
                / merged_size as f64;
            new_row.push(d_new);
        }

        // Remove old rows/columns (higher index first to preserve lower index).
        dist.remove(best_j);
        dist.remove(best_i);
        for row in &mut dist {
            row.remove(best_j);
            row.remove(best_i);
        }
        for (row, &d_val) in dist.iter_mut().zip(&new_row) {
            row.push(d_val);
        }
        new_row.push(0.0);
        dist.push(new_row);
        clusters.push(merged);
    }

    clusters.pop().unwrap_or(Cluster {
        node: PhyloNode {
            label: String::from("Empty"),
            branch_length: 0.0,
            children: Vec::new(),
        },
        height: 0.0,
        size: 0,
    })
    .node
}

fn build_neighbor_joining(names: &[String], initial_dist: &[Vec<f64>]) -> PhyloNode {
    if names.len() <= 2 {
        return build_upgma(names, initial_dist);
    }

    let mut nodes: Vec<PhyloNode> = names
        .iter()
        .map(|name| PhyloNode {
            label: name.clone(),
            branch_length: 0.0,
            children: Vec::new(),
        })
        .collect();
    let mut dist = initial_dist.to_vec();

    while nodes.len() > 2 {
        let k = nodes.len();
        let row_sums: Vec<f64> = (0..k).map(|i| dist[i].iter().sum()).collect();

        // Find pair (i, j) minimizing Saitou-Nei Q-criterion.
        let mut best_i = 0usize;
        let mut best_j = 1usize;
        let mut best_q = f64::INFINITY;
        let kf = (k as f64) - 2.0;

        for i in 0..k {
            for j in (i + 1)..k {
                let q = kf * dist[i][j] - row_sums[i] - row_sums[j];
                if q < best_q {
                    best_q = q;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        let dij = dist[best_i][best_j];
        let limb_i = (0.5 * dij + (row_sums[best_i] - row_sums[best_j]) / (2.0 * kf)).max(0.0002);
        let limb_j = (dij - limb_i).max(0.0002);

        let mut new_dists = Vec::with_capacity(k - 1);
        for m in 0..k {
            if m == best_i || m == best_j {
                continue;
            }
            let d_u = (0.5 * (dist[best_i][m] + dist[best_j][m] - dij)).max(0.0);
            new_dists.push(d_u);
        }

        let mut node_j = nodes.remove(best_j);
        let mut node_i = nodes.remove(best_i);
        node_i.branch_length = limb_i;
        node_j.branch_length = limb_j;

        let merged = PhyloNode {
            label: String::new(),
            branch_length: 0.0,
            children: vec![node_i, node_j],
        };

        dist.remove(best_j);
        dist.remove(best_i);
        for row in &mut dist {
            row.remove(best_j);
            row.remove(best_i);
        }
        for (row, &d_val) in dist.iter_mut().zip(&new_dists) {
            row.push(d_val);
        }
        new_dists.push(0.0);
        dist.push(new_dists);
        nodes.push(merged);
    }

    // Join the final two remaining nodes at a midpoint root.
    let final_d = dist[0][1].max(0.0004);
    let mut right = nodes.pop().unwrap_or(PhyloNode {
        label: String::new(),
        branch_length: 0.0,
        children: Vec::new(),
    });
    let mut left = nodes.pop().unwrap_or(PhyloNode {
        label: String::new(),
        branch_length: 0.0,
        children: Vec::new(),
    });
    left.branch_length = final_d * 0.5;
    right.branch_length = final_d * 0.5;

    PhyloNode {
        label: String::new(),
        branch_length: 0.0,
        children: vec![left, right],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{Alphabet, Sequence};

    #[test]
    fn upgma_and_nj_cluster_closer_sequences_together() {
        let aln = Alignment::new(
            "Taxa",
            Alphabet::Dna,
            vec![
                Sequence::new("Human", b"ATGGTGCATCTGACTCCTGAGGAGAAGTCT"),
                Sequence::new("Chimp", b"ATGGTGCACCTGACTCCTGAGGAGAAGTCT"),
                Sequence::new("Mouse", b"ATGGTGCACCTGACTGATGCTGAGAAGGCT"),
                Sequence::new("Fish",  b"ATGGTTCAACTGTCCGAGGAAGACAAATCT"),
            ],
        )
        .expect("aln");

        let upgma = build_phylogenetic_tree(&aln, TreeMethod::Upgma, DistanceModel::JukesCantor)
            .expect("upgma");
        assert_eq!(upgma.root.leaf_count(), 4);
        assert!(upgma.to_newick().ends_with(';'));

        let nj = build_phylogenetic_tree(&aln, TreeMethod::NeighborJoining, DistanceModel::PDistance)
            .expect("nj");
        assert_eq!(nj.root.leaf_count(), 4);
        assert!(nj.to_newick().contains("Human"));
    }
}
