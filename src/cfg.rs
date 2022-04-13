use crate::proposal;
use crate::tree;
use crate::params;

use crate::Engine;
use crate::util::{PriorDist, log_normal_categories};

use rand::Rng;
use rand::seq::SliceRandom;

use statrs::distribution::{Gamma, ContinuousCDF};

#[derive(Debug)]
pub enum TreePrior {
    Uniform { root: PriorDist },
    Coalescent { num_intervals: usize },
}

impl TreePrior {
    pub fn draw(&self, engine: &mut Engine) -> params::TreePriorParams {
        match self {
            Self::Uniform { root } => {
                params::TreePriorParams::Uniform {
                    root: root.draw(engine),
                }
            },
            Self::Coalescent { num_intervals } => { 
                params::TreePriorParams::Coalescent {
                    // dummy value to be replaced when we know the interval count
                    sizes: (0..*num_intervals).map(|_| 1).collect::<Vec<usize>>(),
                    pops: (0..*num_intervals).map(|_| 1.0).collect::<Vec<f64>>(),
                }
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct Calibration {
    pub low: f64,
    pub high: f64,
}

#[derive(Debug)]
pub struct Constraint {
    pub tips: Vec<usize>,
    pub ancestor: Option<usize>,
}

#[derive(Debug)]
pub struct TreeModel {
    pub prior: TreePrior,
    pub data: tree::TreeData,
    pub calibrations: Vec<(usize, Calibration)>,
    pub constraints: Vec<Constraint>,
}

impl TreeModel {
    pub fn draw(&self, engine: &mut Engine) -> params::TreeParams {
        let mut prior = self.prior.draw(engine);
        let mut root_node = tree::TreeNode::blank();

        match prior {
            params::TreePriorParams::Uniform { root } => {
                root_node.height = root;
                root_node.length = 1.0;
            },
            params::TreePriorParams::Coalescent { .. } => {
                root_node.height = self.calibrations.iter().map(|(_, c)| c.high).fold(0.0, |a: f64, b: f64| a.max(b)) * 1.1;
                root_node.length = 1.0;
            },
        }

        let mut nodes: Vec<tree::TreeNode> = vec![root_node];
        let mut active = vec![]; 
        for tip in &self.data.tips {
            nodes.push(tree::TreeNode { id: tip.id, parent: 0, lchild: 0, rchild: 0,
                                        length: 0.0, height: 0.0, clade: false });
            active.push(tip.id);
        }

        for (id, cal) in &self.calibrations {
            let dist = PriorDist::Uniform { low: cal.low, high: cal.high };
            nodes[*id].height = dist.draw(engine);
        }

        loop {
            if let [l, r] = active[..] {
                nodes[0].lchild = l;
                nodes[0].rchild = r;

                nodes[l].parent = 0;
                nodes[r].parent = 0;

                let llen = nodes[0].height - nodes[l].height;
                let rlen = nodes[0].height - nodes[r].height;

                nodes[l].length = llen;
                nodes[r].length = rlen;
                break
            }

            let mut low = 0;
            let mut low_height = f64::MAX;
            for (i, n) in active.iter().enumerate() {
                let h = nodes[*n].height;
                if h < low_height {
                    low = i;
                    low_height = h; 
                }
            }

            let a = active[low];
            active.swap_remove(low);
            let a_height = low_height;

            let b_id = engine.rng.gen_range(0..active.len());
            let b = active.swap_remove(b_id);
            let b_height = nodes[b].height;

            let mut parent = tree::TreeNode::blank();
            parent.id = nodes.len();

            parent.lchild = a;
            parent.rchild = b;

            nodes[a].parent = parent.id;
            nodes[b].parent = parent.id;

            let p = engine.rng.gen_range(2..10) as f64;
            let height = b_height + (nodes[0].height - b_height) * (1.0 / p);
            parent.height = height;

            nodes[a].length = height - a_height;
            nodes[b].length = height - b_height;

            active.push(parent.id);
            nodes.push(parent);
        }

        let tree = tree::Tree { nodes };

        // update coalescent prior with default interval sizes
        match prior {
            params::TreePriorParams::Coalescent { ref mut sizes, .. } => {
                if let TreePrior::Coalescent { num_intervals } = self.prior {

                    let coalescents = tree.coalescents();
                    let some_coals = coalescents.len() / num_intervals;
                    *sizes = vec![];

                    for i in 0 .. num_intervals - 1 {
                        let ceil = some_coals * (i + 1);
                        sizes.push(ceil - 1);
                    }

                    sizes.push(coalescents.len() - 1);
                } else { unreachable!() }
            },
            _ => {}, 
        }

        params::TreeParams {
            prior,
            tree,
        }

    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
        match self.prior {
            TreePrior::Uniform { .. } => {
                let mut high = 0;
                let mut high_val = 0.0;
                let mut low = 0;
                let mut low_val = f64::MAX;
                for tip in 1..=self.data.num_tips() {
                    let h = params.tree.tree.nodes[tip].height;
                    if h >= high_val {
                        high_val = h;
                        high = tip;
                    }
                    if h < low_val {
                        low_val = h;
                        low = tip;
                    }
                }
                assert!(high != low);
                let mut product = 0.0;
                // we are going to skip the highest and lowest tip
                for tip in 1..=self.data.num_tips() {
                    if (tip != high) && (tip != low) {
                        product += (1.0/params.tree.tree.dist(0, tip)).ln();
                    }
                }

                if let params::TreePriorParams::Uniform { ref mut root } = params.tree.prior {
                    *root = params.tree.tree.nodes[0].height; 
                    product -= root.ln();
                } else { unreachable!() }

                product
            },
            TreePrior::Coalescent { num_intervals } => {
                if let params::TreePriorParams::Coalescent { ref sizes, ref pops } = params.tree.prior {
                    let intervals = params.tree.tree.intervals();
                    let coalescents = params.tree.tree.coalescents();

                    let mut product = 0.0;
                    
                    let mut gi_ptr = 0;
                    let mut gi_final = coalescents[sizes[gi_ptr]].2;

                    let mut switch = false;

                    for (i, interval) in intervals.iter().enumerate() {
                        if switch {
                            gi_ptr += 1;
                            gi_final = coalescents[sizes[gi_ptr]].2;
                            switch = false;
                        }

                        let k_i = interval.lineages as f64;
                        let pop = pops[gi_ptr];

                        let k_i_choose_2 = ((k_i * (k_i - 1.0))) / 2.0;

                        let k_i_choose_2_div_pop = k_i_choose_2 / pop;

                        let term1 = if interval.coalescent {
                            k_i_choose_2_div_pop.ln()
                        }
                        else { 0.0 };

                        let term2 = k_i_choose_2_div_pop * interval.width;

                        product += term1 - term2;

                        if i == gi_final {
                            switch = true;
                        }
                    }

                    // population prior

                    product -= pops[0].ln();
                    for j in 1..num_intervals {
                        product -= pops[j-1].ln();
                        product -= pops[j] / pops[j-1];
                    }

                    product

                } else { unreachable!() }
            },
        }
    }

}

#[derive(Debug)]
pub enum SubstitutionModel {
    BinaryGTR { pi_one: PriorDist },
}

impl SubstitutionModel {
    pub fn draw(&self, engine: &mut Engine) -> params::SubstitutionModelParams {
        match self {
            Self::BinaryGTR { pi_one: dist } => {
                params::SubstitutionModelParams::BinaryGTR { pi_one: dist.draw(engine) }
            },
        }
    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
        let pi_one_param = match params.traits.subst {
            params::SubstitutionModelParams::BinaryGTR { pi_one } => pi_one,
        };
        match self {
            SubstitutionModel::BinaryGTR { ref pi_one } => pi_one.log_density(pi_one_param)
        }
    }
}

#[derive(Debug)]
pub struct ASRV {
    pub enabled: bool,
    pub shape: PriorDist,
    pub ncats: usize,
}

impl ASRV {
    pub fn params_for_shape(&self, shape: f64) -> params::ASRVParams {
        let f_ncats = self.ncats as f64; // we c# now

        let alpha = shape;
        let beta = 1.0 / alpha;

        let gamma = Gamma::new(alpha, beta).unwrap();

        let mut cuts = Vec::with_capacity(self.ncats);
        cuts.push(0.0);
        for i in 1..self.ncats {
            let frac = (i as f64) / f_ncats; 
            // inverse_cdf is equivalent to boost's quantile
            // however it is a generic implementation (binary search)
            // and is probably slow and inaccurate
            // would be better to use something bespoke
            cuts.push(gamma.inverse_cdf(frac));
        }

        let gamma_plus = Gamma::new(alpha + 1.0, beta).unwrap();
        let mut plus_cdf_points = Vec::with_capacity(self.ncats + 1);
        plus_cdf_points.push(0.0);
        for i in 1..self.ncats {
            plus_cdf_points.push(gamma_plus.cdf(cuts[i]));
        }
        plus_cdf_points.push(1.0);

        let mut rates = Vec::with_capacity(self.ncats);
        for i in 0..self.ncats {
            let numerator = plus_cdf_points[i+1] - plus_cdf_points[i];
            // denominator is demonstrably always 1/ncats
            rates.push(numerator * f_ncats);
        }

        params::ASRVParams { shape, rates }
    }

    pub fn draw(&self, engine: &mut Engine, _sites: i32) -> params::ASRVParams {
        if !self.enabled {
            let rates = (0..self.ncats).map(|_| 1.0).collect::<Vec<f64>>();
            return params::ASRVParams { shape: 0.0, rates: rates }
        }
        let shape = self.shape.draw(engine);
        self.params_for_shape(shape)
    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
        if self.enabled {
            self.shape.log_density(params.traits.asrv.shape)
        }
        else {
            0.0
        }
    }
}

#[derive(Debug)]
pub struct ABRV {
    pub enabled: bool,
    pub shape: PriorDist, 
}

impl ABRV {
    pub fn draw(&self, engine: &mut Engine, ntips: usize) -> params::ABRVParams {
        if !self.enabled {
            return params::ABRVParams { shape: 0.0,
                                        rates: vec![],
                                        assignment: vec![] };
        }

        let shape = self.shape.draw(engine);
        let rates = log_normal_categories(shape, 2 * ntips - 2);
        let mut assignment = vec![usize::MAX]; // first value should be invalid

        let mut cat_idxs: Vec<usize> = (0..(2 * ntips - 2)).collect();
        cat_idxs.shuffle(&mut engine.rng);
        assignment.extend(cat_idxs);

        params::ABRVParams { shape, rates, assignment }
    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
        if self.enabled {
            self.shape.log_density(params.traits.abrv.shape)
        }
        else {
            0.0
        }
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

impl TraitsModel {
    pub fn draw(&self, engine: &mut Engine, data_traits: i32, num_tips: usize) -> params::TraitsParams {
        params::TraitsParams {
            num_traits: data_traits, //TODO impl num traits prior
            subst: self.subst.draw(engine),
            base: self.base.draw(engine),
            asrv: self.asrv.draw(engine, data_traits),
            abrv: self.abrv.draw(engine, num_tips),
        }
    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
         let log_prior_base = self.base.log_density(params.traits.base); 
         let log_prior_asrv = self.asrv.log_prior_likelihood(params);
         let log_prior_abrv = self.abrv.log_prior_likelihood(params);
         let log_prior_subst = self.subst.log_prior_likelihood(params);

         log_prior_base + log_prior_asrv + log_prior_abrv + log_prior_subst
    }
}

#[derive(Debug)]
pub struct Configuration {
    pub tree: TreeModel,
    pub traits: TraitsModel,
}

impl Configuration {
    pub fn draw(&self, engine: &mut Engine) -> params::Parameters {
        params::Parameters {
            tree: self.tree.draw(engine),
            traits: self.traits.draw(engine, self.tree.data.traits, self.tree.data.num_tips()),
        }
    }

    pub fn log_prior_likelihood(&self, params: &mut params::Parameters) -> f64 {
        self.tree.log_prior_likelihood(params) + self.traits.log_prior_likelihood(params)
    }

    pub fn get_moves(&self) -> proposal::Propose {
        let mut propose = proposal::Propose::empty();

        propose.add_move("Pi One", proposal::PiOneMove::new(), 1);
        propose.add_move("Base Rate", proposal::BaseRateMove::new(), 1);
        propose.add_move("Tree LOCAL", proposal::TreeLocalMove::new(), 3);

        match self.tree.prior {
             TreePrior::Coalescent { .. } => {
                 propose.add_move("Coal Interval Resize", proposal::CoalescentIntervalResize::new(), 1);
                 propose.add_move("Coal Population Rescale", proposal::CoalescentPopulationRescale::new(), 1);
                 propose.add_move("Coal Population Augment", proposal::CoalescentPopulationAugment::new(), 1);
             },
             TreePrior::Uniform { .. } => {},
        }

        if self.tree.calibrations.len() > 0 {
            propose.add_move("Tree Tip", proposal::TreeTipMove::new(), 2);
        }

        if self.traits.asrv.enabled {
            propose.add_move("ASRV Shape", proposal::ASRVShapeMove::new(), 1);
        }

        if self.traits.abrv.enabled {
            propose.add_move("ABRV Shape", proposal::ABRVShapeMove::new(), 1);
            propose.add_move("ABRV Category", proposal::ABRVCategorySwap::new(), 2);
        }

        propose.lock();
        propose
    }
}

