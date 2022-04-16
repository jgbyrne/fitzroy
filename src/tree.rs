use beagle;

use std::collections::HashMap;

use crate::params;
use crate::{Damage, Engine};

#[derive(Clone, Debug)]
pub struct Interval {
    pub width: f64,
    pub end: f64,
    pub lineages: usize,
    pub coalescent: bool,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Relation {
    Equivalent,
    Subset,
    Superset,
    Disjoint,
    Intersecting,
}

static CHUNK_LEN: usize = 8 * std::mem::size_of::<usize>();

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Clade {
    pub ntaxa: usize,
    flags: Vec<usize>,
}

impl Clade {
    pub fn null() -> Self {
        Clade { ntaxa: 0, flags: vec![] }
    }

    pub fn blank(ntaxa: usize) -> Self {
        let chunks = ntaxa / CHUNK_LEN + if (ntaxa % CHUNK_LEN == 0) {0} else {1};
        Clade { ntaxa, flags: (0..chunks).map(|_| 0).collect::<Vec<usize>>() }
    }

    pub fn with(ntaxa: usize, taxon_ids: &Vec<usize>) -> Self {
        let mut new = Self::blank(ntaxa);
        for id in taxon_ids {
            new.add(*id);
        }
        new
    }

    pub fn add(&mut self, id: usize) {
        let chunk = id / CHUNK_LEN;
        let idx = id % CHUNK_LEN;
        self.flags[chunk] |= 1 << idx;
    }

    pub fn relation(&self, other: &Self) -> Relation {
        assert!(self.ntaxa == other.ntaxa);

        let mut may_equiv = true;
        let mut may_super = true;
        let mut may_sub   = true;
        let mut may_disj  = true;

        for i in 0..self.flags.len() {
            if self.flags[i] != other.flags[i] {
                may_equiv = false;
            }

            let and = self.flags[i] & other.flags[i];
            if and != 0 {
                may_disj = false;
            }

            if and != other.flags[i] {
                may_sub = false;
            }

            if and != self.flags[i] {
                may_super = false;
            }
        }

        let rel = if may_equiv {
            Relation::Equivalent
        }
        else if may_super {
            assert!(!may_sub);
            Relation::Superset
        }
        else if may_sub {
            Relation::Subset
        }
        else if may_disj {
            Relation::Disjoint
        }
        else {
            Relation::Intersecting
        };

        //println!("{:?}\n{:064b}{:064b}\n{:064b}{:064b}\n\n", rel, self.flags[0], self.flags[1], other.flags[0], other.flags[1]);

        rel
    }
/*
        let mut rels = vec![];
        for i in 0..self.flags.len() {
            if self.flags[i] == other.flags[i] {
                rels.push(Relation::Equivalent);
            }
            else {
                let and = self.flags[i] & other.flags[i];
                if and == 0 {
                    rels.push(Relation::Disjoint);
                }
                else if and == self.flags[i] {
                    rels.push(Relation::Superset);
                }
                else if and == other.flags[i] {
                    rels.push(Relation::Subset);
                }
                else {
                    return Relation::Intersecting;
                }
            }
        }
        rels.into_iter().fold(Relation::Equivalent, move |acc, rel| {
            match acc {
                Relation::Equivalent => {
                    rel
                },
                a => {
                    if a == rel { a } else { Relation::Intersecting }
                }
            }
        })
    }
*/
    pub fn disjoint(&self, other: &Self) -> bool {
        assert!(self.ntaxa == other.ntaxa);
        for i in 0..self.flags.len() {
            if self.flags[i] & other.flags[i] != 0 {
                return false;
            }
        }
        return true;
    }

    pub fn union(&self, other: &Self) -> Self {
        assert!(self.ntaxa == other.ntaxa);
        let mut new = Self::blank(self.ntaxa);
        for i in 0..new.flags.len() {
            new.flags[i] |= self.flags[i];
            new.flags[i] |= other.flags[i];
        }
        new
    }

    pub fn has(&self, id: usize) -> bool {
        let chunk = id / CHUNK_LEN;
        let idx = id % CHUNK_LEN;
        (self.flags[chunk] & (1 << idx)) != 0
    }
}

#[derive(Clone, Debug)]
pub struct TreeNode {
    pub id: usize,
    pub parent: usize,
    pub lchild: usize,
    pub rchild: usize,
    pub length: f64,
    pub height: f64,
    pub clade: bool,
}

impl TreeNode {
    pub fn blank() -> Self {
        TreeNode {
            id: 0,
            parent: 0,
            lchild: 0,
            rchild: 0,
            length: 0.0,
            height: 0.0,
            clade: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Tree {
    pub ntaxa: usize,
    pub nodes: Vec<TreeNode>,
    pub level_order: Vec<usize>,
    pub clades: Vec<Clade>,
}

impl Tree {
    pub fn new(nodes: Vec<TreeNode>, ntaxa: usize) -> Self {
        let mut tree = Tree {
            ntaxa,
            nodes,
            level_order: vec![],
            clades: vec![],
        };
        tree.update();
        tree
    }

    pub fn is_leaf(&self, id: usize) -> bool {
        self.nodes[id].lchild == 0 && self.nodes[id].rchild == 0
    }

    /*
    pub fn num_leaves(&self) -> usize {
        let mut ctr = 0;
        for node in &self.nodes {
            if node.lchild == 0 && node.rchild == 0 {
               ctr += 1; 
            }
        }
        ctr
    }
    */

    pub fn dist(&self, a: usize, b: usize) -> f64 {
        self.nodes[a].height - self.nodes[b].height
    }

    pub fn length(&self, id: usize) -> f64 {
        self.nodes[id].length
    }

    pub fn lchild(&self, id: usize) -> usize {
        self.nodes[id].lchild 
    }

    pub fn rchild(&self, id: usize) -> usize {
        self.nodes[id].rchild 
    }

    pub fn node_parent(&self, id: usize) -> Option<usize> {
        if id == 0 {
            None
        }
        else {
            Some(self.nodes[id].parent)
        }
    }

    pub fn intervals(&self) -> Vec<Interval> {
        assert!(self.nodes.len() > 2);

        let mut cur = Interval {width: 0.0, end: self.nodes[0].height, lineages: 2, coalescent: true}; 

        let lc_height = self.nodes[self.nodes[0].lchild].height;
        let rc_height = self.nodes[self.nodes[0].rchild].height;

        let mut fringe = if lc_height > rc_height {
            vec![(rc_height, self.nodes[0].rchild), (lc_height, self.nodes[0].lchild)]
        }
        else {
            vec![(lc_height, self.nodes[0].lchild), (rc_height, self.nodes[0].rchild)]
        };

        let mut intervals = vec![];
        while let Some(nxt) = fringe.pop() {
            cur.width = cur.end - nxt.0;
            let lineages = cur.lineages;
            intervals.insert(0, cur);

            if self.is_leaf(nxt.1) {
                cur = Interval {width: 0.0, end: nxt.0, lineages: lineages - 1, coalescent: false};
            }
            else {
                cur = Interval {width: 0.0, end: nxt.0, lineages: lineages + 1, coalescent: true};
            }

            // if we have reached t=0 then stop building intervals
            if nxt.0 == 0.0 {
                break
            }

            let mut fringe_add = |add_h_n: (f64, usize)| {
                let mut added = false;
                for (i, (h, _)) in fringe.iter().enumerate() {
                    if *h > add_h_n.0 {
                        fringe.insert(i, add_h_n);
                        added = true; break; 
                    }
                }
                if !added {
                    fringe.push(add_h_n);
                }
            };

            let lc = self.nodes[nxt.1].lchild;
            let rc = self.nodes[nxt.1].rchild;
            
            if lc != 0 { 
                fringe_add((self.nodes[lc].height, lc));
            }
            if rc != 0 { 
                fringe_add((self.nodes[rc].height, rc));
            }
        }

        intervals 
    }

    pub fn coalescents(&self) -> Vec<(Interval, usize, usize)> {
        let mut coalescents = vec![];
        let mut non_coal_ctr = 0;
        let mut i = 0;
        for int in self.intervals() {
            non_coal_ctr += 1;
            if int.coalescent {
                coalescents.push((int, non_coal_ctr, i));
                non_coal_ctr = 0
            }
            i += 1;
        }

        coalescents
    }

    pub fn update(&mut self) {
        let mut ordering = vec![];
        let mut cur_lvl = vec![0];
        let mut nxt_lvl = vec![];

        while cur_lvl.len() > 0 {
            for node_id in &cur_lvl {
                let lchild = self.nodes[*node_id].lchild;
                let rchild = self.nodes[*node_id].rchild;
                
                if lchild != 0 { nxt_lvl.push(lchild); }
                if rchild != 0 { nxt_lvl.push(rchild); }
            }
            ordering.extend(cur_lvl);
            cur_lvl = nxt_lvl;
            nxt_lvl = vec![];
        }

        let mut clades = (0..self.nodes.len()).map(|_| Clade::null()).collect::<Vec<Clade>>(); 
        for node_id in ordering.iter().rev() {
            if self.is_leaf(*node_id) {
                let mut clade = Clade::blank(self.ntaxa); 
                clade.add(*node_id);
                clades[*node_id] = clade;
            }
            else {
                let left = self.nodes[*node_id].lchild;
                let right = self.nodes[*node_id].rchild;
                assert!(clades[left].disjoint(&clades[right]));
                clades[*node_id] = clades[left].union(&clades[right]);
            }
            //println!("{:?} {:064b}{:064b}", *node_id, clades[*node_id].flags[0],clades[*node_id].flags[1]);
        }

        self.clades = clades;
        self.level_order = ordering;
    }

    pub fn beagle_operations(&self, engine: &Engine, damage: &Damage) -> Vec<beagle::sys::Operation> {
        let inst = engine.beagle();
        let mut op_vec = vec![];
        for node_id in self.level_order.iter().rev() {
            if !self.is_leaf(*node_id)  && damage.is_marked_partials(*node_id) {
                let node   = &self.nodes[*node_id];

                let b_node  = self.beagle_id(*node_id);
                let b_left  = self.beagle_id(node.lchild);
                let b_right = self.beagle_id(node.rchild);

                op_vec.push(
                    inst.build_operation(b_node, b_left, b_right)
                );
            }
        }
        op_vec
    }

    pub fn beagle_edge_updates(&self, base_rate: f64, abrv: Option<&params::ABRVParams>, damage: &Damage) -> Vec<beagle::MatrixUpdate> {
        // skip root (0)
        let mut updates = vec![];
        for tree_node in 1..self.nodes.len() {
            let get_branch_rate = |abrv: &params::ABRVParams| abrv.rates[abrv.assignment[tree_node]];
            if damage.is_marked_matrix(tree_node) {
                updates.push(beagle::MatrixUpdate {
                    model_id: 0,
                    node_id: self.beagle_id(tree_node),
                    edge_length: self.nodes[tree_node].length * base_rate * abrv.map_or(1.0, get_branch_rate),
                });
            }
        }
        updates
    }

    pub fn beagle_id(&self, node_id: usize) -> i32 {
        if node_id == 0 { (self.nodes.len() - 1) as i32 }
        else { (node_id - 1) as i32 }
    }
}

#[derive(Debug)]
pub struct Tip {
    pub id: usize,
    pub name: String,
    pub data: Vec<f64>,
}

#[derive(Debug)]
pub struct Interior {
    pub id: usize,
    pub name: String,
}

#[derive(Debug)]
pub struct TreeData {
    pub traits: i32,
    pub tips: Vec<Tip>,
    pub interiors: HashMap<usize, Interior>,
}

impl TreeData {
    pub fn from_tips(num_traits: i32, tips: Vec<(String, Vec<f64>)>) -> Self {
        let tips = tips.into_iter()
                       .enumerate()
                       .map(|(i, (name, data))| Tip { id: i+1, name, data})
                       .collect::<Vec<Tip>>();
        Self { 
            traits: num_traits,
            tips: tips,
            interiors: HashMap::new()
        }
    }

    pub fn num_tips(&self) -> usize {
        self.tips.len()
    }

    pub fn tip(&self, tip_id: usize) -> &Tip {
        &self.tips[tip_id - 1]
    }

    pub fn tips(&self) -> &Vec<Tip> {
        &self.tips
    }
}

fn write_node_newick(tree: &Tree, data: &TreeData, node: usize) -> String {
    if tree.is_leaf(node) {
        format!("{}:{:.1}", data.tip(node).name, tree.length(node))
    }
    else {
        format!("({},{}):{:.1}", write_node_newick(tree, data, tree.lchild(node)),
                                 write_node_newick(tree, data, tree.rchild(node)),
                                 tree.length(node))
    }
}

pub fn write_newick(tree: &Tree, data: &TreeData) -> String {
    format!("({},{});", write_node_newick(tree, data, tree.lchild(0)), write_node_newick(tree, data, tree.rchild(0)))
}


