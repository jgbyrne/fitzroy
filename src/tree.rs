use beagle;

use std::collections::HashMap;

use crate::cfg::Calibration;

#[derive(Clone, Debug)]
pub struct TreeNode {
    pub id: usize,
    pub parent: usize,
    pub lchild: usize,
    pub rchild: usize,
    pub length: f64,
    pub height: f64,
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
        }
    }
}

#[derive(Clone, Debug)]
pub struct Tree {
    pub nodes: Vec<TreeNode>,
}

impl Tree {
    pub fn is_leaf(&self, id: usize) -> bool {
        self.nodes[id].lchild == 0 && self.nodes[id].rchild == 0
    }

    pub fn num_leaves(&self) -> usize {
        let mut ctr = 0;
        for node in &self.nodes {
            if node.lchild == 0 && node.rchild == 0 {
               ctr += 1; 
            }
        }
        ctr
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

    pub fn edge_lengths(&self) -> Vec<f64> {
        self.nodes.iter().map(|n| n.length).collect::<Vec<f64>>()
    }

    pub fn level_order(&self) -> Vec<usize> {
        let mut ordered = vec![];
        let mut cur_lvl = vec![0];
        let mut nxt_lvl = vec![];

        while cur_lvl.len() > 0 {
            for node_id in &cur_lvl {
                let lchild = self.nodes[*node_id].lchild;
                let rchild = self.nodes[*node_id].rchild;
                
                if lchild != 0 { nxt_lvl.push(lchild); }
                if rchild != 0 { nxt_lvl.push(rchild); }
            }
            ordered.extend(cur_lvl);
            cur_lvl = nxt_lvl;
            nxt_lvl = vec![];
        }
        ordered
    }

    pub fn beagle_operations(&self) -> Vec<beagle::sys::Operation> {
        let mut op_vec = vec![];
        for node_id in self.level_order().iter().rev() {
            let node = &self.nodes[*node_id];
            if node.lchild != 0 && node.rchild != 0 {
                op_vec.push(
                    beagle::sys::Operation {
                        destinationPartials: *node_id as i32,
                        destinationScaleWrite: -1,
                        destinationScaleRead: -1,
                        child1Partials: node.lchild as i32,
                        child1TransitionMatrix: node.lchild as i32,
                        child2Partials: node.rchild as i32,
                        child2TransitionMatrix: node.rchild as i32,
                    }
                );
            }
        }
        op_vec
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


