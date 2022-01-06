use beagle;

use std::collections::HashMap;

use crate::cfg::Calibration;

#[derive(Debug)]
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

#[derive(Debug)]
pub struct Tree {
    pub nodes: Vec<TreeNode>,
}

impl Tree {
    pub fn num_leaves(&self) -> usize {
        let mut ctr = 0;
        for node in &self.nodes {
            if node.lchild == 0 && node.rchild == 0 {
               ctr += 1; 
            }
        }
        ctr
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
    pub data: Vec<f32>,
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
    pub fn from_tips(tips: Vec<(String, Vec<f32>)>) -> Self {
        let n_traits = match tips.get(0) {
            Some((_, data)) => data.len(),
            None => 0,
        };
        let tips = tips.into_iter()
                       .enumerate()
                       .map(|(i, (name, data))| Tip { id: i+1, name, data})
                       .collect::<Vec<Tip>>();
        Self { 
            traits: n_traits as i32,
            tips: tips,
            interiors: HashMap::new()
        }
    }

    pub fn label_interior(&mut self, int: Interior) {
        self.interiors.insert(int.id, int);
    }
}
