extern crate beagle;
//use rand::Rng;

#[cfg(test)]
mod test;

pub mod cfg;
pub mod params;
pub mod tree;
//pub mod newick;
pub mod nexus;

pub struct Engine {
    inst: beagle::Instance,
}

pub struct MCMC {
    config: cfg::Configuration,
    params: params::Parameters,
    engine: Engine,
}

impl MCMC {
    pub fn new(config: cfg::Configuration) -> Self {

    }
}


