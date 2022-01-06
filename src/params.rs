use crate::tree;
use beagle;

#[derive(Clone)]
pub enum TreePriorParams { 
    Uniform { root: f64 },
}

#[derive(Clone)]
pub struct TreeParams {
    pub prior: TreePriorParams,
    pub tree: tree::Tree,
}

#[derive(Clone)]
pub enum SubstitutionModelParams {
    BinaryGTR { pi_one: f64 },
}

#[derive(Clone)]
pub struct ASRVParams {
    pub shape: f64,
    pub rates: Vec<f64>,
}

#[derive(Clone)]
pub struct ABRVParams {
    pub shape: f64,
}

#[derive(Clone)]
pub struct TraitsParams {
    pub num_traits: i32,
    pub subst: SubstitutionModelParams,
    pub base: f64,
    pub asrv: ASRVParams,
    pub abrv: ABRVParams,
}

impl TraitsParams {
    pub fn model(&self) -> beagle::Model {
        let (freqs, evals, evecs, inv_evecs) = match self.subst {
            SubstitutionModelParams::BinaryGTR { pi_one } => {
                let pi_zero = 1.0 - pi_one;

                let freqs = vec![pi_zero, pi_one];

                // project notebook 1 page 51
                let scale = 1.0 / (2.0 * pi_zero * pi_one);

                let a = pi_one * scale;
                let b = pi_zero * scale;

                let evals = vec![0.0, -a-b];
                let evecs = vec![1.0, 1.0, -(a/b), 1.0];
                let inv_evecs = vec![pi_zero, -pi_zero, pi_one, pi_zero];

                (freqs, evals, evecs, inv_evecs)
            },
        };

        beagle::Model {
            state_freqs: freqs,
            eigenvalues: evals,
            eigenvectors: evecs,
            inv_eigenvectors: inv_evecs,
            category_rates: vec![1.0, 1.0, 1.0, 1.0], // TODO ASRV
            category_probs: vec![0.25, 0.25, 0.25, 0.25],
        }
    }
}

#[derive(Clone)]
pub struct Parameters {
    pub tree: TreeParams,
    pub traits: TraitsParams,
}


