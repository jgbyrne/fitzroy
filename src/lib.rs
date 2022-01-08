extern crate beagle;
use rand::Rng;
use rand::rngs::ThreadRng;

#[cfg(test)]
mod test;

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

    fn instance(&mut self) -> &mut beagle::Instance {
        if !self.run { panic!("Tried to access beagle instance but not running") }
        self.inst.as_mut().unwrap()
    }
}

pub struct MCMC<'c> {
    pub config: &'c cfg::Configuration,
    pub params: params::Parameters,
    pub engine: Engine,
    pub last_log_likelihood: f64, 
}

impl<'c> MCMC<'c> {
    pub fn new(config: &'c cfg::Configuration) -> Self {
        let mut engine = Engine::forge(&config);
        let params = config.draw(&mut engine);
        
        if engine.start(&config, &params) { println!("Engine started") };

        Self {
            config,
            params,
            engine,
            last_log_likelihood: f64::NEG_INFINITY,
        }
    }

    pub fn clone(&self) -> Self {
        let mut engine = Engine::forge(&self.config);
        engine.start(&self.config, &self.params);

        Self {
            config: self.config,
            params: self.params.clone(),
            engine,
            last_log_likelihood: self.last_log_likelihood,
        }
    }

    pub fn step(&mut self) {
    }

    pub fn log_likelihood(&mut self, damage: &Damage) -> f64 {
        self.engine.instance().set_models(&vec![self.params.traits.model()]);

        let updates = self.params.tree.tree.beagle_edge_updates(damage);
        self.engine.instance().update_matrices(updates);

        let ops = self.params.tree.tree.beagle_operations(&mut self.engine, damage);
        self.engine.instance().perform_operations(ops);
        self.engine.instance().calculate_root_log_likelihood(self.params.tree.tree.beagle_id(0), 0)
    }
}


