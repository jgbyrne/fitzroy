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
        let mut inst = beagle::Instance::new(2, n_sites, 4, n_nodes, n_tips, true, vec![model]);

        for tip in &config.tree.data.tips {
            inst.set_tip_data_partial(tip.id as i32, tip.data.clone());
        }

        inst.update_matrices(0, params.tree.tree.edge_lengths());

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
    pub log_likelihood: f64, 
}

impl<'c> MCMC<'c> {
    pub fn new(config: &'c cfg::Configuration) -> Self {
        let mut engine = Engine::forge(&config);
        let params = config.draw(&mut engine);
        
        engine.start(&config, &params);

        Self {
            config,
            params,
            engine,
            log_likelihood: 0.0,
        }
    }

    pub fn clone(&self) -> Self {
        let mut engine = Engine::forge(&self.config);
        engine.start(&self.config, &self.params);

        Self {
            config: self.config,
            params: self.params.clone(),
            engine,
            log_likelihood: self.log_likelihood,
        }
    }

    pub fn log_likelihood(&mut self) -> f64 {
        let ops = self.params.tree.tree.beagle_operations();
        self.engine.instance().perform_operations(ops);
        self.engine.instance().calculate_root_log_likelihood(0, 0)
    }
}


