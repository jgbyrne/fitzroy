extern crate beagle;
use rand::Rng;
use rand::rngs::ThreadRng;

#[cfg(test)]
mod test;

pub mod util;
pub mod proposal;
pub mod cfg;
pub mod params;
pub mod tree;
//pub mod newick;
pub mod nexus;

pub struct Damage {
    partials: Vec<bool>,
    matrices: Vec<bool>,
}

impl Damage {
    pub fn blank(tree: &tree::Tree) -> Self {
        let n_nodes = tree.nodes.len();
        Damage { 
            partials: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
            matrices: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
        }
    }

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

    pub fn mark_partials_to_root(&mut self, tree: &tree::Tree, node_id: usize) {
        self.mark_partials(node_id);
        let mut cur = node_id;
        while let Some(parent) = tree.node_parent(cur) {
            self.mark_partials(parent);
            cur = parent;
        }
    }

    pub fn mark_all_matrices(&mut self) {
        self.matrices.iter_mut().for_each(|item| *item = true);
    }
}

pub struct Engine {
    run: bool,
    rng: ThreadRng,
    inst: Option<beagle::Instance>,
}

impl Engine {
    fn forge<'c>(config: &'c cfg::Configuration) -> Self {
        Engine { run: false, rng: rand::thread_rng(), inst: None }
    }

    fn start<'c>(&mut self, config: &'c cfg::Configuration, params: &params::Parameters) -> bool {
        if self.run { panic!("Tried to start running Engine") }

        let model = params.traits.model();
        let n_sites = config.tree.data.traits;
        let n_tips = config.tree.data.num_tips() as i32;
        let n_nodes = (n_tips * 2) - 1; //TODO trial buffers
        let mut inst = beagle::Instance::new(2, n_sites, 4, n_nodes, n_tips, 1,
                                             true, true, false);

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

pub struct MCMC<'c> {
    pub config: &'c cfg::Configuration,
    pub params: params::Parameters,
    pub propose: proposal::Propose,
    pub engine: Engine,
    pub last_log_likelihood: f64, 
}

impl<'c> MCMC<'c> {
    pub fn new(config: &'c cfg::Configuration) -> Self {
        let mut engine = Engine::forge(&config);
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

    pub fn clone(&self) -> Self {
        let mut engine = Engine::forge(&self.config);
        engine.start(&self.config, &self.params);

        Self {
            config: self.config,
            params: self.params.clone(),
            propose: self.config.get_moves(),
            engine,
            last_log_likelihood: self.last_log_likelihood,
        }
    }

    pub fn step(&mut self) {
        let result = self.propose.make_move(&self.config, &mut self.params, &mut self.engine);

        let accepted = if result.log_prior_likelihood_delta != f64::NEG_INFINITY {
            let likelihood = self.log_likelihood_with_damage(&result.damage);
            let likelihood_delta = likelihood - self.last_log_likelihood;
            let log_ratio = result.log_prior_likelihood_delta + likelihood_delta + result.log_hastings_ratio;

            if log_ratio > util::log_uniform(&mut self.engine) {
                self.last_log_likelihood = likelihood;
                true
            } else { false }
        } else { false };

        if !accepted { (result.revert)(self); }
    }

    pub fn log_likelihood_with_damage(&mut self, damage: &Damage) -> f64 {
        let tree = &self.params.tree.tree;

        let updates = tree.beagle_edge_updates(damage);
        let ops = tree.beagle_operations(&mut self.engine, damage);
        
        let inst = self.engine.beagle_mut();
        inst.set_models(&vec![self.params.traits.model()]);
        inst.update_matrices(updates);
        inst.perform_operations(ops);

        inst.calculate_root_log_likelihood(tree.beagle_id(0), 0)
    }

    pub fn log_likelihood(&mut self) -> f64 {
        self.log_likelihood_with_damage(&Damage::full(&self.params.tree.tree))
    }
}


