use crate::Engine;

use rand::Rng;
use rand::seq::SliceRandom;
use rand::distributions::Distribution;

use statrs::distribution::{Exp, Continuous};

#[derive(Debug)]
pub enum PriorDist {
    Reciprocal,
    Uniform { low: f64, high: f64 },
    Exponential { l: f64 },
}

impl PriorDist {
    pub fn draw(&self, engine: &mut Engine) -> f64 {
        match self {
            PriorDist::Reciprocal => {
                unimplemented!();
            },
            PriorDist::Uniform { low, high } => {
                let r: f64 = engine.rng.gen();
                low + (r * (high - low))
            },
            PriorDist::Exponential { l } => {
                let dist = Exp::new(*l).unwrap();
                dist.sample(&mut engine.rng)
            },
        }
    }

    pub fn log_density(&self, engine: &mut Engine, x: f64) -> f64 {
        match self {
            PriorDist::Reciprocal => {
                unimplemented!();
            },
            PriorDist::Uniform { low, high } => {
                (1.0 / (high - low)).ln()
            },
            PriorDist::Exponential { l } => {
                let dist = Exp::new(*l).unwrap();
                dist.ln_pdf(x)
            },
        }
    }
}

pub fn log_uniform(engine: &mut Engine) -> f64 {
    PriorDist::Uniform { low: 0.0, high: 1.0 }.draw(engine).ln()
}


