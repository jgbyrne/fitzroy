use crate::{MCMC, Engine, Damage};
use crate::params;
use crate::cfg;
use crate::util::{PriorDist, log_normal_categories};

use std::boxed::Box;
use rand::Rng;
use rand::seq::SliceRandom;

pub struct Propose {
    moves: Vec<Box<dyn Move>>,
    weights: Vec<usize>,
    probs: Vec<f32>,
    locked: bool,
}

impl Propose {
    pub fn empty() -> Self {
        Self {
            moves: vec![],
            weights: vec![],
            probs: vec![],
            locked: false,
        }
    }

    pub fn add_move(&mut self, mv: Box<dyn Move>, weight: usize) {
        assert!(!self.locked);
        self.moves.push(mv);
        self.weights.push(weight);
    }

    pub fn lock(&mut self) {
        assert!(!self.locked);
        assert!(self.moves.len() > 0);
        let wsum = self.weights.iter().sum::<usize>() as f32;
        self.probs = self.weights.iter()
                                 .map(|w| (*w as f32) / wsum)
                                 .collect::<Vec<f32>>();
        self.locked = true;
    }

    pub fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        assert!(self.locked);
        let target: f32 = engine.rng.gen();
        let mut acc = 0.0;
        for i in 0..self.moves.len() {
            acc += self.probs[i];
            if acc > target {
                return self.moves[i].make_move(config, params, engine);
            }
        }
        unreachable!();
    }
}

pub type Revert<'c> = dyn FnOnce(&'c mut MCMC) -> ();

pub struct MoveResult<'c> {
    pub log_prior_likelihood_delta: f64,
    pub log_hastings_ratio: f64,
    pub damage: Damage,
    pub revert: Box<Revert<'c>>,
}

pub trait Move {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c>;
}

// Morph the tree topology using the LOCAL move for rooted trees

pub struct TreeLocalMove { }
impl TreeLocalMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for TreeLocalMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        let cur_log_prior_likelihood = config.tree.log_prior_likelihood(params);

        let tree = &mut params.tree.tree;

        // choose an internal edge at random
        assert!(tree.nodes.len() > 3);
        let u = loop {
            let n = engine.rng.gen_range(1..tree.nodes.len());
            let u = tree.node_parent(n).unwrap();
            if u != 0 {
                break u;
            }
        };

        let u_clade = tree.nodes[u].clade;

        let v = tree.node_parent(u).unwrap();

        let a = tree.lchild(u);
        let b = tree.rchild(u);

        let (c, u_on_left) = if tree.lchild(v) == u {
            (tree.rchild(v), true)
        } else {
            (tree.lchild(v), false)
        };

        let a_b_c = vec![a, b, c];

        let nodes_backup = vec![tree.nodes[a].clone(),
                                tree.nodes[b].clone(),
                                tree.nodes[c].clone(),
                                tree.nodes[u].clone(),
                                tree.nodes[v].clone()];

        let mut damage = Damage::blank(&tree);
        damage.mark_partials_to_root(&tree, a);
        damage.mark_partials(b);
        damage.mark_partials(c);
        damage.mark_matrix(a);
        damage.mark_matrix(b);
        damage.mark_matrix(c);
        damage.mark_matrix(u);

        // case: v is not the root
        let hastings = if v != 0 {
            damage.mark_matrix(v);

            let w = tree.node_parent(v).unwrap();
            let a_dist = tree.dist(w, a);
            let b_dist = tree.dist(w, b);
            let c_dist = tree.dist(w, c);
            let u_dist = tree.dist(w, u);
            
            let dists = vec![a_dist, b_dist, c_dist];

            let mut h_ord = vec![(0, a_dist), (1, b_dist), (2, c_dist)];
            h_ord.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            let (h1, h2, h3) = (h_ord[0].0, h_ord[1].0, h_ord[2].0);

            let h1_interval = PriorDist::Uniform { low: 0.0, high: dists[h1] };
            let h2_interval = PriorDist::Uniform { low: 0.0, high: dists[h2] };

            let y = h1_interval.draw(engine);
            let x = h2_interval.draw(engine);

            let (u_star, v_star) = if x > y {
                (x, y)
            }
            else {
                (y, x)
            };

            // u, v heights and edges 
            tree.nodes[v].length = v_star;
            tree.nodes[u].length = u_star - v_star;

            tree.nodes[u].height = tree.nodes[w].height - u_star;
            tree.nodes[v].height = tree.nodes[w].height - v_star;

            if u_star < dists[h1] {
                // any of {a, b, c} can join to v
                let v_child = engine.rng.gen_range(0..3);
                tree.nodes[a_b_c[v_child]].parent = v;
                tree.nodes[a_b_c[v_child]].length = dists[v_child] - v_star;

                if u_on_left {
                    tree.nodes[v].rchild = a_b_c[v_child];
                }
                else {
                    tree.nodes[v].lchild = a_b_c[v_child];
                }

                // the other two join to u
                let u_lchild = (v_child + 1) % 3;
                let u_rchild = (v_child + 2) % 3;

                tree.nodes[a_b_c[u_lchild]].parent = u;
                tree.nodes[a_b_c[u_lchild]].length = dists[u_lchild] - u_star;

                tree.nodes[a_b_c[u_rchild]].parent = u;
                tree.nodes[a_b_c[u_rchild]].length = dists[u_rchild] - u_star;
               
                tree.nodes[u].lchild = a_b_c[u_lchild];
                tree.nodes[u].rchild = a_b_c[u_rchild];

                // hastings ratio
                if u_dist > c_dist { 3.0 } else { 1.0 } 
            }
            else {
                // child at h1 must join to v
                tree.nodes[a_b_c[h1]].parent = v;
                tree.nodes[a_b_c[h1]].length = dists[h1] - v_star;

                if u_on_left {
                    tree.nodes[v].rchild = a_b_c[h1];
                }
                else {
                    tree.nodes[v].lchild = a_b_c[h1]
                }

                tree.nodes[a_b_c[h2]].parent = u;
                tree.nodes[a_b_c[h2]].length = dists[h2] - u_star;

                tree.nodes[a_b_c[h3]].parent = u;
                tree.nodes[a_b_c[h3]].length = dists[h3] - u_star;
               
                tree.nodes[u].lchild = a_b_c[h2];
                tree.nodes[u].rchild = a_b_c[h3];

                // hastings ratio
                if u_dist < c_dist { 1.0 / 3.0 } else { 1.0 }
            }
        }
        else {
            // v is the root
            let a_dist = tree.dist(v, a);
            let b_dist = tree.dist(v, b);
            let c_dist = tree.dist(v, c);
            let u_dist = tree.length(u);

            let dists = vec![a_dist, b_dist, c_dist];

            let mut h_ord = vec![(0, a_dist), (1, b_dist), (2, c_dist)];
            h_ord.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            let (h1, h2, h3) = (h_ord[0].0, h_ord[1].0, h_ord[2].0);

            let base: f64 = 1.5;
            let scale = base.powf(PriorDist::Uniform {low: 0.0, high: 1.0}.draw(engine) - 0.5);

            let mut dists_star = vec![0.0, 0.0, 0.0]; 
            dists_star[h1] = dists[h1] * scale;
            
            let delta = dists_star[h1] - dists[h1];
            dists_star[h2] = dists[h2] + delta;
            dists_star[h3] = dists[h3] + delta;

            let ratio = dists_star[h1] / dists[h1];

            let u_star = PriorDist::Uniform {low: 0.0, high: dists_star[h2]}.draw(engine);

            // u, v heights and edges 
            tree.nodes[v].height += delta;
            tree.nodes[u].length = u_star;
            tree.nodes[u].height = tree.nodes[v].height - u_star;

            if u_star < dists_star[h1] {
                 // any of {a, b, c} can join to v
                let v_child = engine.rng.gen_range(0..3);
                tree.nodes[a_b_c[v_child]].parent = v;
                tree.nodes[a_b_c[v_child]].length = dists_star[v_child];

                if u_on_left {
                    tree.nodes[v].rchild = a_b_c[v_child];
                }
                else {
                    tree.nodes[v].lchild = a_b_c[v_child];
                }

                // the other two join to u
                let u_lchild = (v_child + 1) % 3;
                let u_rchild = (v_child + 2) % 3;

                tree.nodes[a_b_c[u_lchild]].parent = u;
                tree.nodes[a_b_c[u_lchild]].length = dists_star[u_lchild] - u_star;

                tree.nodes[a_b_c[u_rchild]].parent = u;
                tree.nodes[a_b_c[u_rchild]].length = dists_star[u_rchild] - u_star;
               
                tree.nodes[u].lchild = a_b_c[u_lchild];
                tree.nodes[u].rchild = a_b_c[u_rchild];
                
                // hastings ratio
                if u_dist > c_dist { 3.0 * ratio } else { ratio }
            }
            else {
                // child at h1 must join to v
                tree.nodes[a_b_c[h1]].parent = v;
                tree.nodes[a_b_c[h1]].length = dists_star[h1];

                if u_on_left {
                    tree.nodes[v].rchild = a_b_c[h1];
                }
                else {
                    tree.nodes[v].lchild = a_b_c[h1]
                }

                tree.nodes[a_b_c[h2]].parent = u;
                tree.nodes[a_b_c[h2]].length = dists_star[h2] - u_star;

                tree.nodes[a_b_c[h3]].parent = u;
                tree.nodes[a_b_c[h3]].length = dists_star[h3] - u_star;
               
                tree.nodes[u].lchild = a_b_c[h2];
                tree.nodes[u].rchild = a_b_c[h3];

                // hastings ratio
                if u_dist < c_dist { 1.0 / 3.0 } else { 1.0 }
            }
        };

        let revert = move |chain: &mut MCMC| {
            for node in nodes_backup {
                let id = node.id;
                chain.params.tree.tree.nodes[id] = node;
            }
        };

        let log_prior_likelihood_delta = if u_clade && tree.nodes[c].parent != v {
            f64::NEG_INFINITY
        } else { 
            config.tree.log_prior_likelihood(params) - cur_log_prior_likelihood
        };

        MoveResult {
            log_prior_likelihood_delta,
            log_hastings_ratio: hastings,
            damage,
            revert: Box::new(revert),
        }
    }
}


// The LOCAL move for rooted trees doesn't ever adjust tip heights
// This is fine for `present-day` tips, but we need the ability to adjust calibrated tips

pub struct TreeTipMove { }
impl TreeTipMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for TreeTipMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        assert!(config.tree.calibrations.len() > 0); // This move should only be added if we have calibrations
        let (node_id, cfg::Calibration { low, high }) = *config.tree.calibrations.choose(&mut engine.rng).unwrap();

        let cur_log_prior_likelihood = config.tree.log_prior_likelihood(params);

        let cur_height = params.tree.tree.nodes[node_id].height;
        let cur_length = params.tree.tree.nodes[node_id].length;

        let window = (high - low) / 20.0;
        let proposal = PriorDist::Uniform { low: cur_height - window, high: cur_height + window };
        let mut new_height = proposal.draw(engine);

        if new_height > high {
            new_height = low + (new_height - high);
        }
        else if new_height < low {
            new_height = high + (new_height - low);
        }
        
        let parent_id = params.tree.tree.nodes[node_id].parent;
        params.tree.tree.nodes[node_id].height = new_height;
        params.tree.tree.nodes[node_id].length = params.tree.tree.dist(parent_id, node_id);

        let revert = move |chain: &mut MCMC| {
            let node = &mut chain.params.tree.tree.nodes[node_id];
            node.height = cur_height;
            node.length = cur_length;
        };

        let mut damage = Damage::blank(&params.tree.tree);
        damage.mark_partials_to_root(&params.tree.tree, node_id);
        damage.mark_matrix(node_id);

        MoveResult {
            log_prior_likelihood_delta: config.tree.log_prior_likelihood(params) - cur_log_prior_likelihood,
            log_hastings_ratio: 0.0,
            damage,
            revert: Box::new(revert),
        }

    }
}

// ABRV swap two rate categories

pub struct ABRVCategorySwap {}
impl ABRVCategorySwap {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for ABRVCategorySwap {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        assert!(config.traits.abrv.enabled);

        let cur_assignment = params.traits.abrv.assignment.clone();
        let node_max = 2 * config.tree.data.num_tips() - 1;
        
        let node1 = engine.rng.gen_range(1..node_max); // 1..=2n-2 
        let mut node2 = node1;
        while node2 == node1 {
            node2 = engine.rng.gen_range(1..node_max);
        }

        params.traits.abrv.assignment.swap(node1, node2);

        let revert = move |chain: &mut MCMC| {
            chain.params.traits.abrv.assignment = cur_assignment;
        };

        let mut damage = Damage::blank(&params.tree.tree);
        damage.mark_matrix(node1);
        damage.mark_matrix(node2);
        damage.mark_partials_to_root(&params.tree.tree, node1);
        damage.mark_partials_to_root(&params.tree.tree, node2);

        MoveResult {
            log_prior_likelihood_delta: 0.0,
            log_hastings_ratio: 0.0,
            damage,
            revert: Box::new(revert),
        }
    }
}


// Adjust the ASRV shape 

pub struct ASRVShapeMove { }
impl ASRVShapeMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for ASRVShapeMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        assert!(config.traits.asrv.enabled);

        let cur_asrv_shape = params.traits.asrv.shape;

        let lambda = match config.traits.asrv.shape {
            PriorDist::Exponential { l }  => l,
            _ => unimplemented!(),
        };

        let nudge = 1.0 / lambda;
        let proposal = PriorDist::Uniform { low: cur_asrv_shape - nudge,
                                            high: cur_asrv_shape + nudge };

        let new_asrv_shape = proposal.draw(engine);

        let revert = move |chain: &mut MCMC| {
            chain.params.traits.asrv.shape = cur_asrv_shape;
        };

        let log_prior_likelihood_delta = if new_asrv_shape > 0.0 {
            params.traits.asrv = config.traits.asrv.params_for_shape(new_asrv_shape);
            let old = config.traits.asrv.shape.log_density(cur_asrv_shape);
            let new = config.traits.asrv.shape.log_density(new_asrv_shape);
            new - old
        } else { f64::NEG_INFINITY }; // busted prior

        MoveResult {
            log_prior_likelihood_delta,
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}

// Adjust the ABRV shape 

pub struct ABRVShapeMove { }
impl ABRVShapeMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for ABRVShapeMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        assert!(config.traits.abrv.enabled);

        let cur_abrv_shape = params.traits.abrv.shape;
        let cur_abrv_rates = params.traits.abrv.rates.clone();

        let lambda = match config.traits.abrv.shape {
            PriorDist::Exponential { l }  => l,
            _ => unimplemented!(),
        };

        let nudge = 1.0 / lambda;
        let proposal = PriorDist::Uniform { low: cur_abrv_shape - nudge,
                                            high: cur_abrv_shape + nudge };
        let new_abrv_shape = proposal.draw(engine);

        let n_edges = 2 * config.tree.data.num_tips() - 2;
        params.traits.abrv.rates = log_normal_categories(new_abrv_shape, n_edges);

        let revert = move |chain: &mut MCMC| {
            chain.params.traits.abrv.shape = cur_abrv_shape;
            chain.params.traits.abrv.rates = cur_abrv_rates;
        };

        let log_prior_likelihood_delta = if new_abrv_shape > 0.0 {
            let old = config.traits.abrv.shape.log_density(cur_abrv_shape);
            let new = config.traits.abrv.shape.log_density(new_abrv_shape);
            new - old
        } else { f64::NEG_INFINITY }; // busted prior

        MoveResult {
            log_prior_likelihood_delta,
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}



// Adjust the trait evolution base rate

pub struct BaseRateMove { }
impl BaseRateMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for BaseRateMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        let cur_base_rate = params.traits.base;

        let lambda = match config.traits.base {
            PriorDist::Exponential { l }  => l,
            _ => unimplemented!(),
        };

        let nudge = 1.0 / lambda;
        let proposal = PriorDist::Uniform { low: cur_base_rate - nudge,
                                            high: cur_base_rate + nudge };
        let new_base_rate = proposal.draw(engine);

        params.traits.base = new_base_rate;

        let revert = move |chain: &mut MCMC| {
            chain.params.traits.base = cur_base_rate;
        };

        let log_prior_likelihood_delta = if new_base_rate > 0.0 {
            let old = config.traits.base.log_density(cur_base_rate);
            let new = config.traits.base.log_density(new_base_rate);
            new - old
        } else { f64::NEG_INFINITY }; // busted prior

        MoveResult {
            log_prior_likelihood_delta,
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}

// Adjust the substitution model equilibrium frequencies (for binary GTR) 

pub struct PiOneMove { }
impl PiOneMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for PiOneMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        let cur_pi_one = match params.traits.subst {
            params::SubstitutionModelParams::BinaryGTR { pi_one } => pi_one,
        };

        let (high, low) = match config.traits.subst {
            cfg::SubstitutionModel::BinaryGTR { pi_one: PriorDist::Uniform { high, low } } => (high, low),
            _ => unimplemented!(),
        };

        let window = (high - low) / 20.0;
        let proposal = PriorDist::Uniform { low: cur_pi_one - window, high: cur_pi_one + window };
        let mut new_pi_one = proposal.draw(engine);

        if new_pi_one > high {
            new_pi_one = low + (new_pi_one - high);
        }
        else if new_pi_one < low {
            new_pi_one = high + (new_pi_one - low);
        }

        match params.traits.subst {
            params::SubstitutionModelParams::BinaryGTR { ref mut pi_one } => *pi_one = new_pi_one,
        }

        let revert = move |chain: &mut MCMC| {
            match chain.params.traits.subst {
                params::SubstitutionModelParams::BinaryGTR { ref mut pi_one } => *pi_one = cur_pi_one,
            }
        };

        MoveResult {
            log_prior_likelihood_delta: 0.0, // we happen to know we have a uniform prior
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}
