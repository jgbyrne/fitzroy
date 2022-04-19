// =-=-=-=-= lib.rs =-=-=-==-=
// This is the core of fitzroy, a phylolinguistic MCMC engine

extern crate beagle;

use std::collections::{HashMap, VecDeque};
use rand::Rng;
use rand::rngs::ThreadRng;

#[cfg(test)]
mod test;

pub mod util;
pub mod proposal;
pub mod cfg;
pub mod params;
pub mod tree;
pub mod nexus;

// `Summary` stores snapshots of trees
// :: `snapshots` stores the number of snapshots 
// :: `clades` tracks the frequency of clades
pub struct Summary {
    pub snapshots: usize,
    pub clades: HashMap<tree::Clade, (usize, Vec<f64>)>,
    pub trees: Vec<tree::Tree>
}

impl Summary {
    // Create an empty `Summary`
    pub fn blank() -> Self {
        Summary {
            snapshots: 0,
            clades: HashMap::new(),
            trees: Vec::new(),
        }
    }

    // Snapshot a `Tree`, adding it to this `Summary`
    pub fn snapshot(&mut self, tree: &tree::Tree) {
        self.snapshots += 1;
        for (node_id, clade) in tree.clades.iter().enumerate() {
            self.clades.entry(clade.clone())
                .and_modify(|count_ages| {
                    count_ages.0 += 1;
                    count_ages.1.push(tree.nodes[node_id].height);
                }).or_insert((1, vec![tree.nodes[node_id].height]));
        }
        self.trees.push(tree.clone());
    }

    // Build an maximum-clade-credibility (MCC) tree from this `Summary`
    pub fn mcc_tree(&self) -> tree::Tree {
        // Search every snapshotted tree for the one with the highest credibility
        let mut mcc = self.trees.iter().max_by_key(|tree|    
            tree.clades.iter().map(|c| self.clades[c].0).sum::<usize>()
        ).unwrap().clone();

        // Set each clade's height to its median snapshotted value
        for (node_id, clade) in mcc.clades.iter().enumerate() {
            let mut heights = self.clades[clade].1.clone();
            heights.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let n = heights.len();
            let median_height = match n % 2 {
                0 => {
                    (heights[n / 2] + heights[(n / 2) + 1]) / 2.0
                },
                1 => {
                    heights[(n / 2) + 1]
                },
                _ => unreachable!(),
            };
            mcc.nodes[node_id].height = median_height;
        }

        // Calculate branch lengths for new node heights
        // :: technically they could be negative!
        // :: this is just part of how MCC works
        let mut stack = VecDeque::new();
        stack.push_back(mcc.nodes[0].lchild);
        stack.push_back(mcc.nodes[0].rchild);

        while let Some(n) = stack.pop_front() {
            mcc.nodes[n].length = mcc.dist(mcc.nodes[n].parent, n);
            if mcc.nodes[n].lchild != 0 { stack.push_back(mcc.nodes[n].lchild) }
            if mcc.nodes[n].rchild != 0 { stack.push_back(mcc.nodes[n].rchild) }
        }

        mcc
    }

}

// `Damage` is used to track what state an MCMC move has invalidated
pub struct Damage {
    partials: Vec<bool>, // A flag for every node partial
    matrices: Vec<bool>, // A flag for every node matrix
}

impl Damage {
    // Create `Damage` with all flags `false`
    pub fn blank(tree: &tree::Tree) -> Self {
        let n_nodes = tree.nodes.len();
        Damage { 
            partials: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
            matrices: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
        }
    }

    // Create `Damage` with all flags `true`
    pub fn full(tree: &tree::Tree) -> Self {
        let mut damage = Self::blank(tree);
        damage.mark_all_partials();
        damage.mark_all_matrices();
        damage
    }

    pub fn is_marked_partials(&self, node_id: usize) -> bool {
        self.partials[node_id]
    }

    pub fn is_marked_matrix(&self, node_id: usize) -> bool {
        self.matrices[node_id]
    }

    pub fn mark_partials(&mut self, node_id: usize) {
        self.partials[node_id] = true;
    }

    pub fn mark_matrix(&mut self, node_id: usize) {
        self.matrices[node_id] = true;
    }

    pub fn mark_all_partials(&mut self) {
        self.partials.iter_mut().for_each(|item| *item = true);
    }

    pub fn mark_all_matrices(&mut self) {
        self.matrices.iter_mut().for_each(|item| *item = true);
    }

    // If a node's partials are invalidated, all node partials between
    // it and the root are necessarily also invalidated
    pub fn mark_partials_to_root(&mut self, tree: &tree::Tree, node_id: usize) {
        self.mark_partials(node_id);
        let mut cur = node_id;
        while let Some(parent) = tree.node_parent(cur) {
            self.mark_partials(parent);
            cur = parent;
        }
    }
}

// `Engine` is a struct for the 'implementation state' of the MCMC process
// :: `run` is a flag indicating if the `Engine` has been initialised
// :: `rng` is a random number generator
// :: `partial_damage`: whether the beagle `Instance` supports buffer-flipping
// :: `inst` the beagle `Instance` provided by libhmsbeagle-rs
pub struct Engine {
    run: bool,
    rng: ThreadRng,
    partial_damage: bool,
    inst: Option<beagle::Instance>,
}

impl Engine {
    // Create a new `Engine`
    fn forge<'c>(_config: &'c cfg::Configuration, partial_damage: bool) -> Self {
        Engine { run: false, rng: rand::thread_rng(), partial_damage, inst: None }
    }

    // Intialise the `Engine` 
    fn start<'c>(&mut self, config: &'c cfg::Configuration, params: &params::Parameters) -> bool {
        if self.run { panic!("Tried to start running Engine") }

        // Create the beagle `Instance`
        let n_sites = config.tree.data.traits;
        let n_tips = config.tree.data.num_tips() as i32;
        let n_internals = (n_tips * 2) - 1; /* By the FA Cup Theorem :-) */
        let mut inst = beagle::Instance::new(2, n_sites, 4, n_internals, n_tips, 1,
                                             true, true, self.partial_damage);

        // Assign tip partials
        for tip in &config.tree.data.tips {
            let tip_beagle_id = params.tree.tree.beagle_id(tip.id);
            inst.set_tip_data_partial(tip_beagle_id, tip.data.clone());
        }

        self.inst = Some(inst);
        self.run = true;

        true
    }

    fn beagle(&self) -> &beagle::Instance {
        if !self.run { panic!("Tried to access beagle instance but not running") }
        self.inst.as_ref().unwrap()
    }

    fn beagle_mut(&mut self) -> &mut beagle::Instance {
        if !self.run { panic!("Tried to access beagle instance but not running") }
        self.inst.as_mut().unwrap()
    }
}

// `MCMC` - Core structure for an instance of a Markov Chain Monte Carlo
pub struct MCMC<'c> {
    pub config: &'c cfg::Configuration,
    pub params: params::Parameters,
    pub propose: proposal::Propose,
    pub engine: Engine,
    pub last_log_likelihood: f64, 
}

impl<'c> MCMC<'c> {
    // Create a new MCMC instance for a given `Configuration`
    pub fn new(config: &'c cfg::Configuration) -> Self {
        let mut engine = Engine::forge(&config, true);
        let params = config.draw(&mut engine);
        let propose = config.get_moves();
        
        engine.start(&config, &params);

        let mut chain = Self {
            config,
            params,
            propose,
            engine,
            last_log_likelihood: f64::NEG_INFINITY,
        };
        chain.last_log_likelihood = chain.log_likelihood();
        chain
    }

    // Clone a new instance from this one
    // (Intended for multi-chain support, not actually used as of yet)
    pub fn clone(&self) -> Self {
        let mut engine = Engine::forge(&self.config, self.engine.partial_damage);
        engine.start(&self.config, &self.params);

        Self {
            config: self.config,
            params: self.params.clone(),
            propose: self.config.get_moves(),
            engine,
            last_log_likelihood: self.last_log_likelihood,
        }
    }

    // If we have enabled partial damage, we can flip onto new buffers in the beagle `Instance`
    pub fn flip_by_damage(&mut self, damage: &Damage) {
        assert!(self.engine.partial_damage);
        let inst = self.engine.beagle_mut();
        for node in 0..self.params.tree.tree.nodes.len() {
            let beagle_id = self.params.tree.tree.beagle_id(node);
            if damage.is_marked_partials(node) { inst.flip_alt_partials(beagle_id) }
            if node != 0 && damage.is_marked_matrix(node) { inst.flip_alt_matrix(beagle_id); }
        }
    }

    // Perform one step of the Metropolis-Hastings algorithm
    pub fn step(&mut self) -> bool {

        // Draw a move from the proposal distribution and perform it 
        let (result, move_id) = self.propose.make_move(&self.config, &mut self.params, &mut self.engine);

        // If we support partial damage then flip the appropriate buffers 
        if self.engine.partial_damage { self.flip_by_damage(&result.damage); }

        let accepted = if result.log_prior_likelihood_delta != f64::NEG_INFINITY {

            // Calculate likelihood of data given paramaterisation
            let likelihood = if self.engine.partial_damage {
                self.log_likelihood_with_damage(&result.damage)
            } else {
                self.log_likelihood()
            };

            if likelihood > 0.0 { 
                // Something has gone very wrong!
                // :: We could try and revert but since this is probably a fatal bug
                // :: we might as well leave the system state alone to help debugging
                // if self.engine.partial_damage { self.flip_by_damage(&result.damage); }
                // (result.revert)(self);
                return false; 
            }

            // Calculate acceptance probability
            let likelihood_delta = likelihood - self.last_log_likelihood;
            let log_ratio = result.log_prior_likelihood_delta + likelihood_delta + result.log_hastings_ratio;

            // Accept if acceptance probability beats a uniform random variable in [0, 1]
            if log_ratio > util::log_uniform(&mut self.engine) {
                self.last_log_likelihood = likelihood;
                true
            } else { false }

        // If the prior likelihood delta is NEG_INFINITY then reject out-of-hand
        } else { false };

        if !accepted { 
            // Revert to pre-move state 
            if self.engine.partial_damage { self.flip_by_damage(&result.damage); }
            (result.revert)(self);
        }
        else {
            // Tally number of acceptances per move
            self.propose.moves[move_id].3 += 1;
        }

        // MCMC step completed without error
        true 
    }
    
    // Calculate root log likelihood for some damage 
    pub fn log_likelihood_with_damage(&mut self, damage: &Damage) -> f64 {
        let tree = &self.params.tree.tree;

        let abrv = if self.config.traits.abrv.enabled {
            Some(&self.params.traits.abrv)
        } else {
            None
        };

        let updates = tree.beagle_edge_updates(self.params.traits.base, abrv, damage);
        let ops = tree.beagle_operations(&mut self.engine, damage);
       
        let inst = self.engine.beagle_mut();
        inst.set_models(&vec![self.params.traits.model()]);
        inst.update_matrices(updates);
        inst.perform_operations(ops);

        inst.calculate_root_log_likelihood(tree.beagle_id(0), 0)
    }

    // Fully calculate root log likelihood (full damage)
    pub fn log_likelihood(&mut self) -> f64 {
        self.log_likelihood_with_damage(&Damage::full(&self.params.tree.tree))
    }

    // Fully calculate prior likelihood
    pub fn log_prior_likelihood(&mut self) -> f64 {
        self.config.log_prior_likelihood(&mut self.params)
    }

    // Fully calculate posterior likelihood
    pub fn log_posterior_likelihood(&mut self) -> f64 {
        self.log_likelihood() + self.log_prior_likelihood()
    }
}
