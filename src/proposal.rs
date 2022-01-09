use crate::{MCMC, Engine, Damage};
use crate::params;
use crate::cfg;
use crate::util::PriorDist;

use std::boxed::Box;
use rand::Rng;

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

pub type Revert<'c> = dyn Fn(&'c mut MCMC) -> ();

pub struct MoveResult<'c> {
    pub log_prior_likelihood_delta: f64,
    pub log_hastings_ratio: f64,
    pub damage: Damage,
    pub revert: Box<Revert<'c>>,
}

pub trait Move {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c>;
}

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
            log_prior_likelihood_delta: 0.0,
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}

/*
pub struct RootAgeMove { }
impl RootAgeMove {
    pub fn new() -> Box<Self> {
        Box::new(Self { })
    }
}

impl Move for RootAgeMove {
    fn make_move<'c>(&self, config: &cfg::Configuration, params: &mut params::Parameters, engine: &mut Engine) -> MoveResult<'c> {
        let cur_root = match params.tree.prior {
            params::TreePriorParams::Uniform { root } => root,
        };

        let (high, low) = match config.tree.prior {
            cfg::TreePrior::Uniform { root: PriorDist::Uniform { high, low } } => (high, low),
            _ => unimplemented!(),
        };

        let window = (high - low) / 20.0;
        let proposal = PriorDist::Uniform { low: cur_root - window, high: cur_root + window };
        let mut new_root = proposal.draw(engine);

        if new_root > high {
            new_root = low + (new_root - high);
        }
        else if new_root < low {
            new_root = high + (new_root - low);
        }

        println!("{} {}", cur_root, new_root);

        match params.tree.prior {
            params::TreePriorParams::Uniform { ref mut root } => { *root = new_root; },
        };

        let revert = move |chain: &mut MCMC| {
            match chain.params.tree.prior {
                params::TreePriorParams::Uniform { ref mut root } => { *root = cur_root; },
            }
        };

        MoveResult {
            log_prior_likelihood_delta: 0.0,
            log_hastings_ratio: 0.0,
            damage: Damage::full(&params.tree.tree),
            revert: Box::new(revert),
        }
    }
}
*/

