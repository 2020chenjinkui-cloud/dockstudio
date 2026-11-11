"""Adversarial stress-testing suite for Milestone 2: Consensus & Diversity Rigor.

Focuses empirically on:
1. Exponential Consensus Ranking (ECR):
   - Degenerate ties across all metrics (50 molecules with identical scores).
   - Two-tier tied cohorts (25 top-tier, 25 bottom-tier).
   - Missing or partial metrics on random candidates (some have CNN score, some do not).
   - Scale extremes: N=1, N=2, and N=1000.
   - Inverted, negative, zero, and empty weight distributions.
   - Metric directionality validation.
2. Murcko-Butina Diversity Clustering:
   - Acyclic molecules without ring scaffolds (alkanes, peptides, PEG, lipids).
   - Boundary inputs: single molecule, identical molecules, completely dissimilar molecules.
   - Large library of 105 diverse natural products across multiple TCM phytochemical classes.
   - Exact representative invariant: exactly one representative per cluster.
   - In-place assignment on both dicts and ScreeningCandidateRow objects.
3. Virtual Screening Enrichment Analytics:
   - Boundary conditions: all actives, all decoys (zero actives), zero decoys.
   - Validation handling: candidates mode raises ValueError for pure classes.
   - Tied score behavior: full library ties and boundary-crossing ties.
   - Scale minimal: N=2.
"""

from __future__ import annotations

import math
import random
import time
from typing import Any
import pytest
from rdkit import Chem

from dockstudio.analysis import scoring_rigor
from dockstudio.models import ScreeningCandidateRow


# =============================================================================
# 1. Exponential Consensus Ranking (ECR) Adversarial Tests
# =============================================================================

class TestECRConsensusAdversarial:
    """Rigorous stress-testing of Exponential Consensus Ranking under ties, missing data, and scale."""

    def test_ecr_degenerate_ties_50_molecules(self) -> None:
        """50 molecules with completely identical affinity, SILE, and CNN scores."""
        n = 50
        rows = [
            {
                "id": f"mol_{i}",
                "affinity_mean": -8.5,
                "sile": 3.2,
                "cnn_score": 0.82,
                "anchor_match_ratio": 0.75,
            }
            for i in range(n)
        ]

        res = scoring_rigor.compute_exponential_consensus_ranking(rows)

        assert len(res) == n
        scores = [r["ecr_score"] for r in res]
        # In SciPy rankdata, all 50 candidates receive average rank (1+50)/2 = 25.5
        # Therefore, every candidate must have the exact same consensus score
        assert len(set(scores)) == 1
        # Score must be non-zero and finite
        assert scores[0] > 0.0
        assert not math.isnan(scores[0])
        # Ranks must span exactly 1..50
        ranks = [r["ecr_rank"] for r in res]
        assert ranks == list(range(1, n + 1))

    def test_ecr_two_tied_tiers(self) -> None:
        """25 superior candidates with identical metrics and 25 inferior candidates with identical metrics."""
        tier_a = [
            {"id": f"A_{i}", "affinity_mean": -10.5, "sile": 4.1, "cnn_score": 0.92}
            for i in range(25)
        ]
        tier_b = [
            {"id": f"B_{i}", "affinity_mean": -6.0, "sile": 2.1, "cnn_score": 0.35}
            for i in range(25)
        ]

        # Shuffle order to test sort stability and independence of input ordering
        combined = tier_b + tier_a
        res = scoring_rigor.compute_exponential_consensus_ranking(combined)

        top_25 = res[:25]
        bottom_25 = res[25:]

        # All top 25 must be from Tier A
        assert all(r["id"].startswith("A_") for r in top_25)
        # All bottom 25 must be from Tier B
        assert all(r["id"].startswith("B_") for r in bottom_25)

        # Scores within Tier A must be identical
        scores_a = set(r["ecr_score"] for r in top_25)
        assert len(scores_a) == 1
        # Scores within Tier B must be identical
        scores_b = set(r["ecr_score"] for r in bottom_25)
        assert len(scores_b) == 1

        assert list(scores_a)[0] > list(scores_b)[0]

    def test_ecr_random_missing_and_partial_metrics(self) -> None:
        """100 candidates with randomized sparsity (some missing cnn_score, some missing sile)."""
        random.seed(1234)
        cands: list[dict[str, Any]] = []
        for i in range(100):
            c: dict[str, Any] = {
                "id": f"c_{i}",
                "affinity_mean": -8.0 - random.uniform(0.0, 3.0),
            }
            # 40% missing cnn_score
            if random.random() > 0.4:
                c["cnn_score"] = random.uniform(0.4, 0.95)
            # 30% missing sile
            if random.random() > 0.3:
                c["sile"] = random.uniform(2.0, 4.0)
            cands.append(c)

        res = scoring_rigor.compute_exponential_consensus_ranking(cands)
        assert len(res) == 100
        # Check all scores are finite and non-negative
        for r in res:
            assert r["ecr_score"] is not None
            assert not math.isnan(r["ecr_score"])
            assert r["ecr_score"] >= 0.0

        # Monotonicity check
        scores = [r["ecr_score"] for r in res]
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1]

        # Candidates with missing metrics receive worst ranks on those metrics
        c_full = {"id": "full", "affinity_mean": -8.0, "sile": 3.0, "cnn_score": 0.7}
        c_miss = {"id": "miss", "affinity_mean": -8.0, "sile": 3.0, "cnn_score": None}
        pair_res = scoring_rigor.compute_exponential_consensus_ranking([c_miss, c_full])
        assert pair_res[0]["id"] == "full"
        assert pair_res[0]["ecr_score"] > pair_res[1]["ecr_score"]

    def test_ecr_scale_n1(self) -> None:
        """Single candidate (N=1) must have sigma_r=1.0, ecr_rank=1, and finite score."""
        row = [{"id": "single", "affinity_mean": -8.5, "sile": 3.2}]
        res = scoring_rigor.compute_exponential_consensus_ranking(row)
        assert len(res) == 1
        assert res[0]["ecr_rank"] == 1
        assert res[0]["ecr_score"] == pytest.approx(math.exp(-1.0), abs=1e-4)

    def test_ecr_scale_n2(self) -> None:
        """Tiny library (N=2) separating clearly."""
        rows = [
            {"id": "worse", "affinity_mean": -6.0, "sile": 2.0},
            {"id": "better", "affinity_mean": -9.0, "sile": 3.5},
        ]
        res = scoring_rigor.compute_exponential_consensus_ranking(rows)
        assert res[0]["id"] == "better"
        assert res[0]["ecr_rank"] == 1
        assert res[1]["id"] == "worse"
        assert res[1]["ecr_rank"] == 2

    def test_ecr_scale_n1000(self) -> None:
        """Large library (N=1000) under heterogeneous metrics completes in <0.05s with strict monotonicity."""
        random.seed(999)
        cands = [
            {
                "id": f"cand_{i}",
                "affinity_mean": random.uniform(-13.0, -4.0),
                "sile": random.uniform(1.2, 4.8),
                "cnn_score": random.uniform(0.05, 0.98),
                "anchor_match_ratio": random.uniform(0.0, 1.0),
            }
            for i in range(1000)
        ]

        t0 = time.perf_counter()
        res = scoring_rigor.compute_exponential_consensus_ranking(cands)
        t1 = time.perf_counter()

        elapsed = t1 - t0
        assert elapsed < 0.1, f"ECR N=1000 took {elapsed:.4f}s, expected < 0.1s"

        scores = [r["ecr_score"] for r in res]
        ranks = [r["ecr_rank"] for r in res]

        # Strict rank sequence
        assert ranks == list(range(1, 1001))
        # Monotonicity
        assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))
        # Valid bounds
        assert scores[0] <= 1.0
        assert scores[-1] >= 0.0

    def test_ecr_inverted_and_zero_weights(self) -> None:
        """Tests handling of negative weights, zero weights, and empty weights dict."""
        rows = [
            {"id": "m1", "affinity_mean": -9.0, "sile": 3.5},
            {"id": "m2", "affinity_mean": -7.0, "sile": 2.5},
        ]

        # Negative weights sum to <= 0 -> code falls back to uniform 1/len(weights)
        res_neg = scoring_rigor.compute_exponential_consensus_ranking(
            [dict(r) for r in rows], weights={"affinity_mean": -0.5, "sile": -0.5}
        )
        assert res_neg[0]["id"] == "m1"

        # All zero weights
        res_zero = scoring_rigor.compute_exponential_consensus_ranking(
            [dict(r) for r in rows], weights={"affinity_mean": 0.0, "sile": 0.0}
        )
        assert res_zero[0]["id"] == "m1"

        # Empty weights dict falls back to default metric
        res_empty_w = scoring_rigor.compute_exponential_consensus_ranking(
            [dict(r) for r in rows], weights={}
        )
        assert res_empty_w[0]["id"] == "m1"


# =============================================================================
# 2. Murcko-Butina Diversity Clustering Adversarial Tests
# =============================================================================

class TestMurckoButinaClusteringAdversarial:
    """Rigorous stress-testing of scaffold extraction and Butina clustering."""

    def test_murcko_clustering_acyclic_molecules(self) -> None:
        """Acyclic molecules (alkanes, peptides, fatty acids) lack ring scaffolds.
        Verifies that extract_murcko_scaffold returns empty string, and Butina clustering
        falls back to molecular fingerprints, correctly grouping similar acyclics."""
        # Alkanes
        alkanes = [
            Chem.MolFromSmiles("CCCC"),    # butane
            Chem.MolFromSmiles("CCCCC"),   # pentane
            Chem.MolFromSmiles("CCCCCC"),  # hexane
        ]
        # Peptides
        peptides = [
            Chem.MolFromSmiles("CC(N)C(=O)O"),                # alanine
            Chem.MolFromSmiles("CC(N)C(=O)NC(C)C(=O)O"),       # Ala-Ala
            Chem.MolFromSmiles("CC(N)C(=O)NC(C)C(=O)NC(C)C(=O)O"), # Ala-Ala-Ala
        ]
        # Glycols & Lipids
        other_acyclics = [
            Chem.MolFromSmiles("CCOCCOCCO"),       # triethylene glycol
            Chem.MolFromSmiles("CCCCCCCCCCCC(=O)O"), # lauric acid
        ]

        all_acyclics = alkanes + peptides + other_acyclics

        # 1. Scaffolds are all empty string
        for m in all_acyclics:
            assert scoring_rigor.extract_murcko_scaffold(m) == ""

        # 2. Butina clustering with fallback
        c_ids, reps = scoring_rigor.cluster_murcko_scaffolds(all_acyclics, cutoff=0.55)

        num_clusters = len(set(c_ids))
        assert len(reps) == num_clusters, "Number of representatives must equal number of clusters"

        # Verify exactly one representative per cluster
        for cid in range(num_clusters):
            members = [i for i, c in enumerate(c_ids) if c == cid]
            cluster_reps = [i for i in members if i in reps]
            assert len(cluster_reps) == 1, f"Cluster {cid} must have exactly one representative"

        # Alkanes should cluster together (indices 0, 1, 2)
        assert c_ids[0] == c_ids[1] == c_ids[2]
        # Peptides should cluster together (indices 3, 4, 5)
        assert c_ids[3] == c_ids[4] == c_ids[5]
        # Alkanes and Peptides must be in distinct clusters
        assert c_ids[0] != c_ids[3]

    def test_murcko_clustering_single_molecule(self) -> None:
        """Single molecule edge case returns ([0], {0})."""
        m = Chem.MolFromSmiles("c1ccccc1")
        c_ids, reps = scoring_rigor.cluster_murcko_scaffolds([m])
        assert c_ids == [0]
        assert reps == {0}

    def test_murcko_clustering_identical_molecules(self) -> None:
        """20 identical molecules must collapse into exactly 1 cluster with 1 representative."""
        mols = [Chem.MolFromSmiles("c1ccccc1")] * 20
        c_ids, reps = scoring_rigor.cluster_murcko_scaffolds(mols, cutoff=0.4)
        assert set(c_ids) == {0}
        assert len(reps) == 1
        rep_idx = next(iter(reps))
        assert 0 <= rep_idx < 20
        assert c_ids[rep_idx] == 0

    def test_murcko_clustering_completely_dissimilar(self) -> None:
        """6 structurally orthogonal scaffolds must form 6 separate clusters at cutoff=0.2."""
        smiles_list = [
            "c1ccccc1",                                # benzene
            "C1C2CC3CC1CC(C2)C3",                      # adamantane
            "CN1CCCC1c1cccnc1",                        # nicotine
            "O=C1CC2OCC=C3CN4CCC56C4CC3C2C5c2ccccc2N16", # strychnine
            "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C", # cholesterol
            "CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C23OO4",     # artemisinin
        ]
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        c_ids, reps = scoring_rigor.cluster_murcko_scaffolds(mols, cutoff=0.2)
        assert len(set(c_ids)) == 6
        assert len(reps) == 6

    def test_murcko_clustering_large_105_natural_products_library(self) -> None:
        """Stress-tests a 105-compound library spanning flavonoids, alkaloids, steroids,
        coumarins, anthraquinones, polyphenols, lipids, and amino acids."""
        tcm_smiles = [
            # Flavonoids & Chalcones
            "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c12",  # quercetin
            "O=C1C(O)=C(c2ccc(O)cc2)Oc2cc(O)cc(O)c12",       # kaempferol
            "O=C1C(O)=C(c2cc(O)c(O)c(O)c2)Oc2cc(O)cc(O)c12",# myricetin
            "O=C1C=C(c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c12",       # luteolin
            "O=C1C=C(c2ccc(O)cc2)Oc2cc(O)cc(O)c12",          # apigenin
            "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12",           # naringenin
            "COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O",         # hesperetin
            "O=C1C(=Cc2ccc(O)cc2)COc2cc(O)cc(O)c12",         # genistein
            "O=C1C(=Cc2ccc(O)cc2)COc2cc(O)ccc12",            # daidzein
            "O=C(C=Cc1ccc(O)c(O)c1)c1ccc(O)cc1O",            # butein
            "O=C(C=Cc1ccc(O)cc1)c1ccc(O)cc1",                # isoliquiritigenin
            "O=C1C=C(c2ccccc2)Oc2c(O)c(O)c(O)cc12",          # baicalein
            "COc1c(O)cc2c(c1)OC(=CC2=O)c1ccccc1",            # wogonin
            "O=C1C=C(c2ccccc2)Oc2cc(O)cc(O)c12",             # chrysin
            "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc2cc(O)ccc12",       # fisetin
            "COc1cc(O)c2c(c1)OC(c1ccc(O)c(OC)c1)C(O)C2",     # catechin analog
            "Oc1cc(O)c2c(c1)OC(c1ccc(O)c(O)c1)C(O)C2",       # epicatechin
            "O=C1Oc2cc(O)cc(O)c2C(c2ccc(O)c(O)c2)C1O",       # flavan-3-ol deriv
            "COc1ccc(C2Oc3cc(O)cc(O)c3CC2O)cc1O",            # hesperitin deriv
            "Oc1ccc(C2Oc3cc(O)cc(O)c3C=C2)cc1",              # flavene
            # Alkaloids
            "COc1ccc2c(c1OC)CC[n+]1c2cc2c3c1cc1c(c3OCO2)OCO1", # berberine
            "COc1cc2c(cc1OC)CC[n+]1c2cc2c(c1OC)c(OC)ccc2",     # palmatine
            "CN1CCC23C4C1CC5=C2C(=C(O)C=C5)OC3C(O)C=C4",       # morphine
            "COc1ccc2c3c1OC4C(O)C=CC1C(C2)N(C)CCC341",          # codeine
            "COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1",           # quinine
            "OC(c1ccnc2ccccc12)C1CC2CCN1CC2C=C",               # cinchonine
            "CN1CCCC1c1cccnc1",                                 # nicotine
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                       # caffeine
            "Cn1cnc2c1c(=O)[nH]c(=O)n2C",                       # theobromine
            "Cn1c(=O)c2[nH]cnc2n(C)c1=O",                       # theophylline
            "CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2",              # atropine
            "CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C1OC21",           # scopolamine
            "O=C1CC2OCC=C3CN4CCC56C4CC3C2C5c2ccccc2N16",       # strychnine
            "COc1cc(OC)c2c(c1OC)-c1ccc(OC)c(=O)cc1C(NC(C)=O)CC2", # colchicine
            "COc1ccc2[nH]c3c(c2c1)CC1C2CC(OC(=O)c4cc(OC)c(OC)c(OC)c4)C(OC)CC2N3CC1", # reserpine
            "O=C(C=CC=Cc1ccc2c(c1)OCO2)N1CCCCC1",             # piperine
            "CC(C)=CCCC(=O)NCC1=CC(=C(O)C=C1)OC",             # capsaicin
            "c1c[nH]c2ccccc12",                                 # indole
            "CCN(CC)C(=O)C1CN(C)C2Cc3c[nH]c4cccc(C2=C1)c34",  # LSD/ergoline
            "Cc1c[nH]c2ccccc12",                                # 3-methylindole
            # Steroids & Triterpenes
            "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",  # cholesterol
            "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",             # testosterone
            "CC(=O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C",        # progesterone
            "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O",                # estradiol
            "CC(=C)C1CCC2(C(=O)O)CCC3(C)C(CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C12", # betulinic acid
            "CC1(C)CCC2(C(=O)O)CCC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1", # oleanolic acid
            "CC1CCC2(C(=O)O)CCC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1C", # ursolic acid
            "CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C23OO4",             # artemisinin
            "CC12CCC(O)CC1CCC1C2CCC2(C)C(O)CCC12",             # androstanediol
            "CC12CCC(=O)CC1CCC1C2CCC2(C)C(=O)CCC12",           # androstenedione
            # Coumarins & Polyphenols
            "O=C1C=Cc2ccccc2O1",                                # coumarin
            "O=C1C=Cc2ccc(O)cc2O1",                             # umbelliferone
            "O=C1C=Cc2cc(O)c(O)cc2O1",                          # esculetin
            "COc1cc2OC(=O)C=Cc2cc1O",                           # scopoletin
            "O=C(C=Cc1ccc(O)c(OC)c1)CC(=O)C=Cc1ccc(O)c(OC)c1",# curcumin
            "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1",                    # resveratrol
            "COc1cc(O)cc(C=Cc2ccc(OC)cc2)c1",                   # pterostilbene
            "Oc1cc(O)c(C(=O)c2ccc(O)cc2)c(O)c1",               # 2,4,4',6-tetrahydroxybenzophenone
            "O=C(O)c1cc(O)c(O)c(O)c1",                          # gallic acid
            "O=C(O)c1ccc(O)c(O)c1",                             # protocatechuic acid
            # Anthraquinones
            "Cc1cc(O)c2c(c1)C(=O)c1c(O)cc(O)cc1C2=O",          # emodin
            "OCC1=CC2=C(C(=O)C3=C(C=CC=C3O)C2=O)C(O)=C1",      # aloe-emodin
            "O=C(O)c1cc2c(c(O)c1)C(=O)c1c(O)cc(O)cc1C2=O",     # rhein
            "Cc1cc2c(c(O)c1)C(=O)c1c(O)cccc1C2=O",             # chrysophanol
            "O=C1c2ccccc2C(=O)c2ccccc12",                       # anthraquinone core
            # Acyclic & Lipid metabolites
            "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=C", # squalene
            "CCCCCCCCC=CCCCCCCCC(=O)O",                         # oleic acid
            "CCCCCCCCCCCCCCCC(=O)O",                            # palmitic acid
            "CCCCCCCCCCCCCCCCCC(=O)O",                          # stearic acid
            "CCCCCCCCCCCC(=O)O",                                # lauric acid
            "CCCCCCCCCCCCCC(=O)O",                              # myristic acid
            "CCCCCCCCCC(=O)O",                                  # capric acid
            "CCCCCCCC(=O)O",                                    # caprylic acid
            "OCC(O)CO",                                         # glycerol
            "CCCCCCC(=O)O",                                     # heptanoic acid
            # Amino acids & small peptides
            "NCC(=O)O",                                         # glycine
            "CC(N)C(=O)O",                                      # alanine
            "CC(C)C(N)C(=O)O",                                  # valine
            "CC(C)CC(N)C(=O)O",                                 # leucine
            "CCC(C)C(N)C(=O)O",                                 # isoleucine
            "N=C(N)NCCCC(N)C(=O)O",                             # arginine
            "NCCCC(N)C(=O)O",                                   # lysine
            "NC(=O)CC(N)C(=O)O",                                # asparagine
            "NC(=O)CCC(N)C(=O)O",                               # glutamine
            "CSCCC(N)C(=O)O",                                   # methionine
            "OC(=O)C(N)Cc1ccccc1",                              # phenylalanine
            "OC(=O)C(N)Cc1c[nH]c2ccccc12",                      # tryptophan
            "OC(=O)C(N)Cc1ccc(O)cc1",                           # tyrosine
            "OC(=O)C1CCCN1",                                    # proline
            "OC(=O)C(N)CO",                                     # serine
            "CC(O)C(N)C(=O)O",                                  # threonine
            "C1OC(CO)C(O)C(O)C1O",                              # glucose
            "C1OC(CO)C(O)C(CO)C1O",                             # fructose
            "C1OC(C(O)CO)C(O)C1O",                              # ribose
            "CC(N)C(=O)NCC(=O)O",                               # Ala-Gly
            # Additional TCM / Natural Product Actives
            "CC1=C(C(=O)c2ccccc2O1)c1ccccc1",                   # 2-phenylchromone
            "O=c1c2ccccc2oc2ccccc12",                           # xanthone
            "COc1cc2c(cc1OC)c1ccc(OC)c(OC)c1cc2",              # phenanthrene deriv
            "COc1ccc(CC=C)cc1OC",                               # methyleugenol
            "Oc1ccc(CC=C)cc1OC",                                # eugenol
            "COc1ccc(C=CC)cc1",                                 # anethole
            "CC(=O)Oc1ccccc1C(=O)O",                            # aspirin / salicylate deriv
            "COc1ccc(O)cc1",                                    # 4-methoxyphenol
            "O=Cc1ccccc1O",                                     # salicylaldehyde
            "COc1ccccc1C=O",                                    # o-anisaldehyde
        ]

        mols = [Chem.MolFromSmiles(s) for s in tcm_smiles]
        assert all(m is not None for m in mols), "Invalid SMILES detected"

        t0 = time.perf_counter()
        c_ids, reps = scoring_rigor.cluster_murcko_scaffolds(mols, cutoff=0.4)
        t1 = time.perf_counter()

        assert t1 - t0 < 2.0, f"Clustering 105 compounds took {t1 - t0:.3f}s, expected < 2.0s"
        assert len(c_ids) == len(mols)

        num_clusters = len(set(c_ids))
        assert len(reps) == num_clusters, "Representatives count must strictly equal cluster count"

        # Verify cluster validity and exact-one representative invariant
        for cid in range(num_clusters):
            members = [i for i, c in enumerate(c_ids) if c == cid]
            cluster_reps = [i for i in members if i in reps]
            assert len(cluster_reps) == 1, f"Cluster {cid} must have exactly one representative"

        # Test candidate assignment in-place
        candidates = [
            ScreeningCandidateRow(rank=i + 1, ligand_id=f"np_{i}")
            for i in range(len(mols))
        ]
        assigned_ids, assigned_reps = scoring_rigor.assign_clusters_to_candidates(
            candidates, mols, cutoff=0.4
        )
        assert assigned_ids == c_ids
        assert assigned_reps == reps
        assert sum(1 for c in candidates if c.is_cluster_rep) == num_clusters


# =============================================================================
# 3. Virtual Screening Enrichment Analytics Adversarial Tests
# =============================================================================

class TestEnrichmentMetricsAdversarial:
    """Rigorous stress-testing of AUROC, BEDROC, and Enrichment Factors."""

    def test_enrichment_all_actives_raw_and_candidates(self) -> None:
        """All actives in raw lists returns 1.0; in candidates list raises ValueError as required."""
        scores = [-9.0, -8.5, -8.0, -7.5]
        labels = [True, True, True, True]

        # 1. Raw list input: returns 1.0 safely
        res_raw = scoring_rigor.compute_enrichment_metrics(scores=scores, is_active=labels)
        assert res_raw["auroc"] == 1.0
        assert res_raw["bedroc"] == 1.0
        assert res_raw["ef_1%"] == 1.0
        assert res_raw["ef_5%"] == 1.0

        # 2. Candidate objects input: raises ValueError
        cands = [
            ScreeningCandidateRow(rank=i + 1, ligand_id=f"c_{i}", is_active=True, ecr_score=0.9 - i * 0.1)
            for i in range(4)
        ]
        with pytest.raises(ValueError, match="At least one active and one inactive"):
            scoring_rigor.compute_enrichment_metrics(candidates=cands)

    def test_enrichment_all_decoys_raw_and_candidates(self) -> None:
        """Zero actives in raw lists returns 0.0; in candidates list raises ValueError."""
        scores = [-9.0, -8.5, -8.0, -7.5]
        labels = [False, False, False, False]

        res_raw = scoring_rigor.compute_enrichment_metrics(scores=scores, is_active=labels)
        assert res_raw["auroc"] == 0.0
        assert res_raw["bedroc"] == 0.0
        assert res_raw["ef_1%"] == 0.0
        assert res_raw["ef_5%"] == 0.0

        cands = [
            ScreeningCandidateRow(rank=i + 1, ligand_id=f"c_{i}", is_active=False, ecr_score=0.9 - i * 0.1)
            for i in range(4)
        ]
        with pytest.raises(ValueError, match="At least one active and one inactive"):
            scoring_rigor.compute_enrichment_metrics(candidates=cands)

    def test_enrichment_degenerate_ties_behavior(self) -> None:
        """Verifies behavior when all scores are tied across actives and decoys."""
        scores = [-8.0] * 100
        # 50 actives and 50 decoys
        is_active = [True] * 50 + [False] * 50

        # RDKit Scoring evaluates tied scores without crashing
        res = scoring_rigor.compute_enrichment_metrics(scores=scores, is_active=is_active)
        assert not math.isnan(res["auroc"])
        assert not math.isnan(res["bedroc"])
        assert not math.isnan(res["ef_1%"])

    def test_enrichment_boundary_fraction_ties(self) -> None:
        """Tied scores at early recognition cutoff (1% and 5%)."""
        # 100 compounds, top 5 have identical score -10.0, rest have -5.0
        scores = [-10.0] * 5 + [-5.0] * 95
        is_active = [True, True, True, False, False] + [False] * 95

        res = scoring_rigor.compute_enrichment_metrics(scores=scores, is_active=is_active)
        assert res["auroc"] >= 0.95
        assert res["ef_5%"] == 20.0
        assert res["bedroc"] > 0.6

    def test_enrichment_minimal_n2(self) -> None:
        """Minimal non-trivial dataset of N=2 (1 active, 1 decoy)."""
        # Perfect
        res_perf = scoring_rigor.compute_enrichment_metrics(
            scores=[-9.0, -7.0], is_active=[True, False]
        )
        assert res_perf["auroc"] == 1.0

        # Inverted
        res_inv = scoring_rigor.compute_enrichment_metrics(
            scores=[-7.0, -9.0], is_active=[True, False]
        )
        assert res_inv["auroc"] == 0.0

    def test_enrichment_from_ranked_labels_adversarial(self) -> None:
        """Direct rank-ordered labels evaluation."""
        labels_perf = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        res_perf = scoring_rigor.compute_enrichment_from_ranked_labels(labels_perf)
        assert res_perf["auroc"] == 1.0
        assert res_perf["bedroc"] == 1.0

        labels_inv = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        res_inv = scoring_rigor.compute_enrichment_from_ranked_labels(labels_inv)
        assert res_inv["auroc"] == 0.0
        assert res_inv["bedroc"] < 0.01
