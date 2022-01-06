use crate::tree;
use crate::params;
use rand::Rng;

#[derive(Debug)]
pub enum PriorDist {
    Reciprocal,
    Uniform { low: f64, high: f64 },
    Exponential { l: f64 },
}

impl PriorDist {
    pub fn draw(&self) -> f64 {
        match self {
            PriorDist::Reciprocal => {
                unimplemented!();
            },
            PriorDist::Uniform { low, high } => {
                let mut rng = rand::thread_rng();
                let r: f64 = rng.gen();
                low + (r * (high - low))
            },
            PriorDist::Exponential { l } => {
                unimplemented!();
            },
        }
    }
}

#[derive(Debug)]
pub enum TreePrior {
    Uniform { root: PriorDist },
}

impl TreePrior {
    pub fn draw(&self) -> params::TreePriorParams {
        match self {
            Self::Uniform { root } => {
                params::TreePriorParams::Uniform {
                    root: root.draw(),
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Calibration {
    pub low: f64,
    pub high: f64,
}

#[derive(Debug)]
pub struct TreeModel {
    pub prior: TreePrior,
    pub data: tree::TreeData,
    pub calibrations: Vec<(usize, Calibration)>,
}

impl TreeModel {
    pub fn draw(&self) -> params::TreeParams {
        unimplemented!();
    }
}

#[derive(Debug)]
pub enum SubstitutionModel {
    BinaryGTR { pi_one: PriorDist },
}

impl SubstitutionModel {
    pub fn draw(&self) -> params::SubstitutionModelParams {
        match self {
            Self::BinaryGTR { pi_one: dist } => {
                params::SubstitutionModelParams::BinaryGTR { pi_one: dist.draw() }
            },
        }
    }
}

#[derive(Debug)]
pub struct ASRV {
    pub enabled: bool,
    pub shape: PriorDist,
}

impl ASRV {
    pub fn draw(&self, sites: i32) -> params::ASRVParams {
        if !self.enabled {
            return params::ASRVParams { shape: 0.0, rates: vec![] }
        }
        //TODO
        /*let shape = self.shape.draw(); */
        params::ASRVParams { shape: 0.0, rates: vec![] }
    }
}

#[derive(Debug)]
pub struct ABRV {
    pub enabled: bool,
    pub shape: PriorDist, 
}

impl ABRV {
    pub fn draw(&self) -> params::ABRVParams {
        //TODO
        params::ABRVParams { shape: 0.0 }
    }
}

#[derive(Debug)]
pub struct TraitsModel {
    pub num_traits: PriorDist,
    pub subst: SubstitutionModel,
    pub base: PriorDist,
    pub asrv: ASRV,
    pub abrv: ABRV,
}

#[derive(Debug)]
pub struct Configuration {
    pub tree: TreeModel,
    pub traits: TraitsModel,
}
