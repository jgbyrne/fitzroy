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
//    inst: beagle::Instance,
}

pub struct MCMC {
    pub config: cfg::Configuration,
    pub params: params::Parameters,
    pub engine: Engine,
}

impl MCMC {
    pub fn new(config: cfg::Configuration) -> Self {
        let num_traits = config.tree.data.traits;
        let traits = params::TraitsParams {
            num_traits,
            subst: config.traits.subst.draw(),
            base: config.traits.base.draw(),
            asrv: config.traits.asrv.draw(num_traits),
            abrv: config.traits.abrv.draw(),
        };

        let tree = config.tree.draw();

        let params = params::Parameters {
            tree,
            traits,
        };

        Self {
            config,
            params,
            engine: Engine {},
        }
    }
}


