use crate::Engine;

use rand::Rng;
use rand::seq::SliceRandom;
use rand::distributions::Distribution;

use statrs::distribution::{Exp, Continuous};
use statrs::function::erf;

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

    pub fn log_density(&self, x: f64) -> f64 {
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

fn log_normal_quantile(mu: f64, sigma: f64, p: f64) -> f64 {
    let inv_err = erf::erf_inv(2.0 * p - 1.0);
    let exponent = mu + sigma * (std::f64::consts::SQRT_2 * inv_err);
    exponent.exp()
}

pub fn log_normal_categories(sigma: f64, n: usize) -> Vec<f64> {
    let mut cats = Vec::with_capacity(n);

    let mu = -sigma.powi(2) / 2.0;

    let f_n = n as f64;
    for i in 1..=n {
        let f_i = i as f64;
        let p = (f_i - 0.5) / f_n;
        cats.push(log_normal_quantile(mu, sigma, p));
    }

    cats
}
