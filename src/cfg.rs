use crate::proposal;
use crate::tree;
use crate::params;

use crate::Engine;
use crate::util::PriorDist;

use std::boxed::Box;

use rand::Rng;
use rand::seq::SliceRandom;
use rand::distributions::Distribution;

use statrs::distribution::{Exp, Continuous};

#[derive(Debug)]
pub enum TreePrior {
    Uniform { root: PriorDist },
}

impl TreePrior {
    pub fn draw(&self, engine: &mut Engine) -> params::TreePriorParams {
        match self {
            Self::Uniform { root } => {
                params::TreePriorParams::Uniform {
                    root: root.draw(engine),
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Calibration {
    pub low: f64,
    pub high: f64,
}

#[derive(Debug)]
pub struct TreeModel {
    pub prior: TreePrior,
    pub data: tree::TreeData,
    pub calibrations: Vec<(usize, Calibration)>,
}

impl TreeModel {
    pub fn draw(&self, engine: &mut Engine) -> params::TreeParams {
        let prior = self.prior.draw(engine);
        let mut root_node = tree::TreeNode::blank();

        // unsure what to do for priors not conditions on tmrca...
        match prior {
            params::TreePriorParams::Uniform { root } => {
                root_node.height = root;
                root_node.length = 1.0;
            }
        }

        let mut nodes: Vec<tree::TreeNode> = vec![root_node];
        let mut active = vec![]; 
        for tip in &self.data.tips {
            nodes.push(tree::TreeNode { id: tip.id, parent: 0, lchild: 0, rchild: 0,
                                        length: 0.0, height: 0.0});
            active.push(tip.id);
        }

        for (id, cal) in &self.calibrations {
            let dist = PriorDist::Uniform { low: cal.low, high: cal.high };
            nodes[*id].height = dist.draw(engine);
        }

        loop {
            if let [l, r] = active[..] {
                nodes[0].lchild = l;
                nodes[0].rchild = r;

                nodes[l].parent = 0;
                nodes[r].parent = 0;

                let llen = nodes[0].height - nodes[l].height;
                let rlen = nodes[0].height - nodes[r].height;

                nodes[l].length = llen;
                nodes[r].length = rlen;
                break
            }

            let mut low = 0;
            let mut low_height = f64::MAX;
            for (i, n) in active.iter().enumerate() {
                let h = nodes[*n].height;
                if h < low_height {
                    low = i;
                    low_height = h; 
                }
            }

            let a = active[low];
            active.swap_remove(low);
            let a_height = low_height;

            let b_id = engine.rng.gen_range(0..active.len());
            let b = active.swap_remove(b_id);
            let b_height = nodes[b].height;

            let mut parent = tree::TreeNode::blank();
            parent.id = nodes.len();

            parent.lchild = a;
            parent.rchild = b;

            nodes[a].parent = parent.id;
            nodes[b].parent = parent.id;

            let p = engine.rng.gen_range(2..10) as f64;
            let height = b_height + (nodes[0].height - b_height) * (1.0 / p);
            parent.height = height;

            nodes[a].length = height - a_height;
            nodes[b].length = height - b_height;

            active.push(parent.id);
            nodes.push(parent);
        }

        let tree = tree::Tree { nodes };

        params::TreeParams {
            prior,
            tree,
        }

    }
}

#[derive(Debug)]
pub enum SubstitutionModel {
    BinaryGTR { pi_one: PriorDist },
}

impl SubstitutionModel {
    pub fn draw(&self, engine: &mut Engine) -> params::SubstitutionModelParams {
        match self {
            Self::BinaryGTR { pi_one: dist } => {
                params::SubstitutionModelParams::BinaryGTR { pi_one: dist.draw(engine) }
            },
        }
    }
}

#[derive(Debug)]
pub struct ASRV {
    pub enabled: bool,
    pub shape: PriorDist,
}

impl ASRV {
    pub fn draw(&self, engine: &mut Engine, sites: i32) -> params::ASRVParams {
        if !self.enabled {
            return params::ASRVParams { shape: 0.0, rates: vec![] }
        }
        //TODO
        /*let shape = self.shape.draw(); */
        params::ASRVParams { shape: 0.0, rates: vec![] }
    }
}

#[derive(Debug)]
pub struct ABRV {
    pub enabled: bool,
    pub shape: PriorDist, 
}

impl ABRV {
    pub fn draw(&self, engine: &mut Engine) -> params::ABRVParams {
        //TODO
        params::ABRVParams { shape: 0.0 }
    }
}

#[derive(Debug)]
pub struct TraitsModel {
    pub num_traits: PriorDist,
    pub subst: SubstitutionModel,
    pub base: PriorDist,
    pub asrv: ASRV,
    pub abrv: ABRV,
}

impl TraitsModel {
    pub fn draw(&self, engine: &mut Engine, data_traits: i32) -> params::TraitsParams {
        params::TraitsParams {
            num_traits: data_traits, //TODO impl num traits prior
            subst: self.subst.draw(engine),
            base: self.base.draw(engine),
            asrv: self.asrv.draw(engine, data_traits),
            abrv: self.abrv.draw(engine),
        }
    }
}

#[derive(Debug)]
pub struct Configuration {
    pub tree: TreeModel,
    pub traits: TraitsModel,
}

impl Configuration {
    pub fn draw(&self, engine: &mut Engine) -> params::Parameters {
        params::Parameters {
            tree: self.tree.draw(engine),
            traits: self.traits.draw(engine, self.tree.data.traits),
        }
    }

    pub fn get_moves(&self) -> proposal::Propose {
        let mut propose = proposal::Propose::empty();

        propose.add_move(proposal::PiOneMove::new(), 1);
        propose.add_move(proposal::BaseRateMove::new(), 1);

        propose.lock();
        propose
    }
}

