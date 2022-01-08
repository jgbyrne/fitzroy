use super::*;

#[test]
fn test_beagle_partial_tips() {

    let human_str: &str = "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA";
    let chimp_str: &str = "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA";
    let gorilla_str: &str = "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA";

    let partial_seq = |seq: &str| {
        let mut partials: Vec<f64> = vec![];
        for c in seq.chars() {
            match c {
                'A' => partials.extend_from_slice(&[1.0,0.0]),
                'C' => partials.extend_from_slice(&[0.0,1.0]),
                'G' => partials.extend_from_slice(&[0.0,0.0]),
                'T' => partials.extend_from_slice(&[0.0,0.0]),
                _   => partials.extend_from_slice(&[1.0,1.0]),
            }
        }
        partials
    };

    let model = beagle::Model {
        state_freqs: vec![0.5, 0.5],

        eigenvalues: vec![0.0, -2.0],

        eigenvectors: vec![1.0, -1.0, 1.0, 1.0],

        inv_eigenvectors: vec![0.5, 0.5, -0.5, 0.5],

        category_rates: vec![1.0],
        category_probs: vec![1.0],
    };

    let mut inst = beagle::Instance::new(2, 1 as i32, 4, 5, 3, 1, true, false, true);
    inst.set_tip_data_partial(0, vec![0.0, 1.0]);//partial_seq(human_str));
    inst.set_tip_data_partial(1, vec![0.0, 1.0]);//partial_seq(chimp_str));
    inst.set_tip_data_partial(2, vec![1.0, 0.0]);//partial_seq(gorilla_str));

    inst.set_models(&vec![model]);

    let mut updates = vec![];
    let edge_lengths = vec![10000.0, 10000.0, 10000.0, 10000.0];

    for i in 0..edge_lengths.len() {
        updates.push(
            beagle::MatrixUpdate {
                model_id: 0,
                node_id: i as i32,
                edge_length: edge_lengths[i],
            }
        );
    }

    inst.update_matrices(updates);

    let ops = vec![
    beagle::sys::Operation {
        destinationPartials: 3,
        destinationScaleWrite: -1, //OP_NONE
        destinationScaleRead: -1, //OP_NONE
        child1Partials: 0,
        child1TransitionMatrix: 0,
        child2Partials: 1,
        child2TransitionMatrix: 1,
    },
    beagle::sys::Operation {
        destinationPartials: 4,
        destinationScaleWrite: -1, //OP_NONE
        destinationScaleRead: -1, //OP_NONE
        child1Partials: 2,
        child1TransitionMatrix: 2,
        child2Partials: 3,
        child2TransitionMatrix: 3,
    }
    ];

    inst.perform_operations(ops);
    //inst.wait_for_partials(vec![3, 4]);
    println!("Root Sum Log Likelihood: {}", inst.calculate_root_log_likelihood(0, 0));
    println!("Root Sum Log Likelihood: {}", inst.calculate_root_log_likelihood(1, 0));
    println!("Root Sum Log Likelihood: {}", inst.calculate_root_log_likelihood(2, 0));
    println!("Root Sum Log Likelihood: {}", inst.calculate_root_log_likelihood(3, 0));
    println!("Root Sum Log Likelihood: {}", inst.calculate_root_log_likelihood(4, 0));
    assert_eq!(beagle::sys::ReturnCode::SUCCESS, inst.teardown());
}
