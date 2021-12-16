use beagle;

#[derive(Debug)]
pub struct TreeNode {
    pub id: usize,
    pub parent: usize,
    pub lchild: usize,
    pub rchild: usize,
    pub length: f64,
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

pub type Labels = Vec<String>;
