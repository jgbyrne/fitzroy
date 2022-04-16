use crate::Engine;

use rand::Rng;
use rand::distributions::Distribution;

use statrs::distribution::{Exp, Continuous, Gamma, Normal, LogNormal};
use statrs::function::erf;

#[derive(Debug)]
pub enum PriorDist {
    Reciprocal,
    Uniform { low: f64, high: f64 },
    Exponential { l: f64 },
    Gamma { alpha: f64, beta: f64 },
    Normal { mean: f64, sigma: f64 },
    LogNormal { location: f64, scale: f64 },
    HalfNormal { sigma: f64 },
}

impl PriorDist {
    pub fn draw(&self, engine: &mut Engine) -> f64 {
        match self {
            PriorDist::Reciprocal => {
                // chosen by fair dice roll
                // guaranteed to be random
                0.0001
            },
            PriorDist::Uniform { low, high } => {
                let r: f64 = engine.rng.gen();
                low + (r * (high - low))
            },
            PriorDist::Exponential { l } => {
                let dist = Exp::new(*l).unwrap();
                dist.sample(&mut engine.rng)
            },
            PriorDist::Gamma { alpha, beta } => {
                let dist = Gamma::new(*alpha, *beta).unwrap();
                dist.sample(&mut engine.rng)
            },
            PriorDist::Normal { mean, sigma } => {
                let dist = Normal::new(*mean, *sigma).unwrap();
                dist.sample(&mut engine.rng)
            },
            PriorDist::LogNormal { location, scale } => {
                let dist = LogNormal::new(*location, *scale).unwrap();
                dist.sample(&mut engine.rng)
            },
            PriorDist::HalfNormal { sigma } => {
                let dist = Normal::new(0.0, *sigma).unwrap();
                dist.sample(&mut engine.rng).abs()
            },

        }
    }

    pub fn log_density(&self, x: f64) -> f64 {
        match self {
            PriorDist::Reciprocal => {
                - x.ln()
            },
            PriorDist::Uniform { low, high } => {
                (1.0 / (high - low)).ln()
            },
            PriorDist::Exponential { l } => {
                let dist = Exp::new(*l).unwrap();
                dist.ln_pdf(x)
            },
            PriorDist::Gamma { alpha, beta } => {
                let dist = Gamma::new(*alpha, *beta).unwrap();
                dist.ln_pdf(x)
            },
            PriorDist::Normal { mean, sigma } => {
                let dist = Normal::new(*mean, *sigma).unwrap();
                dist.ln_pdf(x)
            },
            PriorDist::LogNormal { location, scale } => {
                let dist = LogNormal::new(*location, *scale).unwrap();
                dist.ln_pdf(x)
            },
            PriorDist::HalfNormal { sigma } => {
                let dist = Normal::new(0.0, *sigma).unwrap();
                dist.ln_pdf(x) + (2.0 as f64).ln()
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
