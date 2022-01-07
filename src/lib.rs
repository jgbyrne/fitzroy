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
    fn blank(tree: &tree::Tree) -> Self {
        let n_nodes = tree.nodes.len();
        Damage { 
            partials: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
            matrices: (0..n_nodes).map(|_| false).collect::<Vec<bool>>(),
        }
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
        let mut inst = beagle::Instance::new(2, n_sites, 4, n_nodes, n_tips,
                                             true, true, false, vec![model]);

        for tip in &config.tree.data.tips {
            let tip_beagle_id = params.tree.tree.beagle_id(tip.id);
            inst.set_tip_data_partial(tip_beagle_id, tip.data.clone());
        }

        let mut edge_lengths = params.tree.tree.beagle_edge_lengths();
        // alternates 
        // edge_lengths.append(&mut edge_lengths.clone());

        inst.update_matrices(0, edge_lengths);

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

    pub fn log_likelihood(&mut self) -> f64 {
        self.params.tree.tree.calculate(&mut self.engine);
        self.engine.instance().calculate_root_log_likelihood(self.params.tree.tree.beagle_id(0), 0)
    }
}


