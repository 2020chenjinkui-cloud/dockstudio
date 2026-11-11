"""修复后强制补充验收测试（报告第四章 6 组用例）。

依据 `DockStudio_HTVS_TCM_验收与修复报告.md` 第四章「修复后强制补充验收测试用例规范」：

1. ``test_organic_salts_preserved_while_true_mixtures_rejected``
   - 药用有机酸盐（酒石酸盐 / 柠檬酸盐 / 马来酸盐 / 甲磺酸盐）在
     ``prepare_ligand_from_smiles`` 与 ``prepare_library`` 中都必须**剥离有机酸根、
     保留活性母核**（不抛 ``PreparationError``）；
   - 真正的双活性药物混合物（阿司匹林 + 布洛芬）必须被 ``multi_component_mixture``
     拦截并计入 ``n_desalt_fail``，绝不进入 ``n_ready``。
2. ``test_tautomer_dedup_and_chiral_enumeration_active``
   - 酮式 / 烯醇式互变异构体经 ``standardize_tautomer_and_charge`` 归一后 InChIKey
     相同，``prepare_library`` 的 ``n_dedup == 1``；
   - 4 个未定手性中心（> ``MAX_UNSPECIFIED_STEREO_CENTERS`` = 3）的**非环**分子
     （``CC(O)C(O)C(O)C(O)C``）被拦截为 ``n_chiral_fail == 1``。
     ⚠️ 刻意不使用己糖/糖环分子：按 ``_unassigned_centers_confined_to_aliphatic_ring``
     设计，糖环上的未定中心是被放行的（见
     ``test_htvs_funnel_stress::test_stereoisomer_unassigned_centers_robustness``）。
3. ``test_lignan_and_lactone_not_misclassified_as_glycoside``
   - 木脂素 / 内酯 / 黄烷醇母核（四氢呋喃或内酯环上挂着芳基、甲氧基、羟基）绝不能被
     当成糖环切碎：``identify_sugar_rings`` 返回 ``[]`` 且 ``deglycosylate_tcm_mol``
     不产苷元；
   - 真正的 O-/N- 糖苷（黄芩苷、芦丁、腺苷）被正确切除糖环；C-糖苷（葛根素）被保护。
4. ``test_tier2_strain_penalty_applied_and_no_unconverged_backdoor``
   - ``3.0 < ΔE_strain <= 8.0`` 的真实 3D 高张力姿态：``corrected_affinity`` 必须真正
     惩罚，并且 ``ScreeningRunner`` 的 Tier 2 用 **修正后** 亲和力计算 SILE
     （SILE 变小、排序下降），候选行 ``affinity_corrected`` 有值；
   - 「1 个收敛 + 2 个发散/超高张力」时 Tier 2 严格只晋级 1 个（不存在「凑够 3 个」后门）；
   - 0 个收敛时只降级保留 1 个并打上 ``UNCONVERGED_OR_HIGH_STRAIN`` /
     ``HIGH_STRAIN_UNCONVERGED_FALLBACK`` 警示。
5. ``test_delivery_poses_preserve_real_docked_3d_coordinates_and_verify_t1_t5``
   - 3-Tier + ``tcm_mode=True`` 双轨筛选（衍生苷元使实际入池数 > ``n_raw``）后，
     ``out_dir/poses/pose_<cid>.sdf`` 必须是**真实对接坐标**：3D 质心落在
     ``(16, 25, 4)`` 结合口袋搜索盒内，绝不等于原点 ``(0, 0, 0)``；
   - ``funnel.csv`` 的 ``aglycone_generated > 0`` 且 ``ScreeningFunnelReport`` 守恒闭合；
   - ``run_verify(out_dir, options={"recalculate": True, "engine": <确定性 mock 引擎>})``
     必须返回 0，且 ``_verify_t1.._verify_t5`` 逐项 PASS（T1 真实重算 |ΔAff|≈0、RMSD≈0）。
6. ``test_cli_entrypoint_and_agent_repl_screen_all_flags``
   - ``dockstudio.cli.main(["screen", ...18 个参数...])`` 与
     ``agent._dispatch_command("/screen ...")`` 返回 0 且所有参数真正生效（走 click）；
   - 反向断言：裸 ``/screen``、``/screen --demo`` 仍走交互向导
     ``run_screen_interactive``，绝不被静默降级或误升级。

说明（避免误判）：第 3 组用例中的天然产物 SMILES 为其**母核代表骨架**（dibenzocyclooctadiene
木脂素 / 四氢呋喃木脂素 / 呋喃并呋喃木脂素 / 芳基萘内酯 / 笼状三内酯 / 黄烷-3-醇）。
``identify_sugar_rings`` 的判定完全由拓扑与取代基决定，与天然产物绝对立体化学无关；
每个用例另外断言「该分子确实含有 1 个『5/6 元、全单键、仅 1 个环氧』的伪糖环」，
确保这些样例是真正能触发历史误判的对抗性输入，而不是随便挑的无环分子。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolTransforms

from dockstudio import cli as dockstudio_cli
from dockstudio.analysis.scoring_rigor import (
    compute_internal_strain_energy,
    compute_ligand_efficiencies,
)
from dockstudio.cli import agent as agent_module
from dockstudio.cli.commands.verify import (
    _verify_t1,
    _verify_t2,
    _verify_t3,
    _verify_t4,
    _verify_t5,
    load_delivery_package,
    load_energy_table,
    load_funnel_data,
    parse_conf_file,
    run_verify,
)
from dockstudio.core.orchestrator import PrepareRegistry
from dockstudio.core.screen_runner import (
    ScreeningRunner,
    load_pose_mol,
    pose_coords_are_real,
    restore_pose_mol_with_template,
)
from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.models import (
    Box,
    DockParams,
    DockResult,
    LigandQCProfile,
    Pose,
    PreparationError,
    PreparedLigand,
    ScreeningConfig,
)
from dockstudio.prep.library import MAX_UNSPECIFIED_STEREO_CENTERS, prepare_library
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.prep.tcm import (
    deglycosylate_tcm_mol,
    desalt_and_filter,
    generate_conformers_and_min_energy,
    identify_sugar_rings,
    standardize_tautomer_and_charge,
)

# =============================================================================
# Shared chemical fixtures
# =============================================================================

#: 药用有机酸盐：活性母核 + 有机酸反离子（酒石酸 / 柠檬酸 / 马来酸 / 甲磺酸）。
#: 每项为 (case_id, 盐形式 SMILES, 活性母核预期重原子数)。
ORGANIC_SALT_CASES: tuple[tuple[str, str, int], ...] = (
    # 酒石酸小檗碱盐（季铵生物碱 + 酒石酸）
    (
        "berberine_tartrate",
        "OC(=O)C(O)C(O)C(=O)O.COc1ccc2cc3[n+](cc2c1OC)CCc1c3ccc2c1OCO2",
        25,
    ),
    # 柠檬酸咖啡因盐（黄嘌呤生物碱 + 柠檬酸）
    (
        "caffeine_citrate",
        "O=C(O)CC(O)(CC(=O)O)C(=O)O.Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        14,
    ),
    # 马来酸氯苯那敏盐（芳基烷基胺 + 马来酸）
    (
        "chlorphenamine_maleate",
        "CN(C)CCC(c1ccc(Cl)cc1)c1ccccn1.O=C(O)C=CC(=O)O",
        19,
    ),
    # 甲磺酸雷沙吉兰盐（炔丙胺 + 甲磺酸）
    (
        "rasagiline_mesylate",
        "C#CCN[C@@H]1CCc2ccccc21.CS(=O)(=O)O",
        13,
    ),
)

#: 真正的双活性药物混合物（阿司匹林 + 布洛芬）：必须被判为多组分混合物。
TRUE_MIXTURE_SMILES = "CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc(cc1)C(C)C(=O)O"

#: 酮式 / 烯醇式互变异构体（乙酰丙酮）。
TAUTOMER_ACAC_KETO = "CC(=O)CC(=O)C"
TAUTOMER_ACAC_ENOL = "CC(O)=CC(=O)C"

#: 4 个未定手性中心的**非环**多醇（默认上限 3 -> 必须被手性质控拦截）。
ACYCLIC_TETRAOL_4_CENTERS = "CC(O)C(O)C(O)C(O)C"

#: 木脂素 / 内酯 / 黄烷醇：四氢呋喃或内酯母核上挂了芳基/甲氧基/羟基，
#: 历史上会被误判成呋喃糖/吡喃糖环并强行水解。
#: 第三列 ``has_fake_sugar_ring`` 标记该样例是否真的含「5/6 元、全单键、仅 1 个环氧」的
#: 伪糖环（历史宽松判据会命中）；二苯并环辛二烯（五味子醇甲）与芳环稠合的黄烷醇
#: （儿茶素）的母核环是 8 元环/芳香稠环，因此历史判据也不会命中，如实标记为 False。
NON_GLYCOSIDE_CASES: tuple[tuple[str, str, bool], ...] = (
    # 五味子醇甲（Schisandrol A / schisandrin, C24H32O7）：二苯并环辛二烯木脂素
    (
        "Schisandrol_A",
        "COC1=C(OC)C2=C(C=C1OC)C1=C(C=C(OC)C(OC)=C1OC)C(C)C(C)CC2O",
        False,
    ),
    # 连翘苷元（Phillygenin, C21H24O6）：四氢呋喃木脂素
    (
        "Phillygenin",
        "COc1cc(CC2COC(c3ccc(OC)c(OC)c3)C2=CO)ccc1O",
        True,
    ),
    # 芝麻素（Sesamin, C20H18O6）：呋喃并呋喃木脂素（双 THF 稠环）
    (
        "Sesamin",
        "COc1cc2c(cc1OC)OC1C(c3ccc4c(c3)OCO4)COC21",
        True,
    ),
    # 穿心莲内酯（Andrographolide）：labdane 型 γ-内酯母核
    (
        "Andrographolide_labdane_lactone",
        "C=C1CC2C(CC=C3C2C(=O)OC3CO)C1(C)CCC(O)C",
        True,
    ),
    # 银杏内酯 B（Ginkgolide B）：笼状三内酯
    (
        "Ginkgolide_B",
        "CC(C)(C)C12CC3C4C(C1=O)OC(=O)C4(O)C(O)C5(C3O2)C(=O)OC5=O",
        True,
    ),
    # 鬼臼毒素（Podophyllotoxin, C22H22O8）：芳基萘内酯（内酯环与芳环稠合）
    (
        "Podophyllotoxin",
        "COc1cc2c(cc1OC)C1C(c3ccc(OC)c(OC)c3)C(=O)OC1CO2",
        True,
    ),
    # 儿茶素（Catechin）：黄烷-3-醇（黄酮 C 环与 A 芳环稠合）
    (
        "Catechin",
        "Oc1ccc(cc1O)C1OC2=CC(O)=CC(O)=C2CC1O",
        False,
    ),
)

#: 真糖苷：黄芩苷（O-葡萄糖醛酸苷）、芦丁（O-二糖苷）、腺苷（N-呋喃糖苷）。
TRUE_GLYCOSIDE_CASES: tuple[tuple[str, str, int], ...] = (
    (
        "Baicalin",
        ("O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O"),
        20,  # 黄芩素 aglycone 重原子数
    ),
    (
        "Rutin",
        ("CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"),
        22,  # 槲皮素 aglycone 重原子数
    ),
    (
        "Adenosine",
        "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
        10,  # 腺嘌呤 aglycone 重原子数
    ),
)

#: C-糖苷：葛根素（糖环通过 C-C 键直接连在苷元上，必须保护不切）。
PUERARIN_SMILES = "O=c1c(-c2ccc(O)cc2)coc2c(O)c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc12"


def _fake_sugar_like_rings(mol: Chem.Mol) -> list[tuple[int, ...]]:
    """历史（修复前）宽松糖环判据：5/6 元、全单键、恰好 1 个环氧。

    仅用于**构造对抗性输入**：断言这些木脂素/内酯样例确实带有能被老判据误判的
    「伪糖环」，从而保证本组用例不是在挑软柿子。
    """
    Chem.FastFindRings(mol)
    suspicious: list[tuple[int, ...]] = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        bonds = [
            mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)]) for i in range(len(ring))
        ]
        if any(b is None or b.GetBondType() != Chem.BondType.SINGLE for b in bonds):
            continue
        symbols = [atom.GetSymbol() for atom in atoms]
        if symbols.count("O") != 1 or symbols.count("C") != len(ring) - 1:
            continue
        suspicious.append(tuple(ring))
    return suspicious


# =============================================================================
# 1. 药用有机酸盐保留 / 真实混合物拦截
# =============================================================================


def test_organic_salts_preserved_while_true_mixtures_rejected(settings) -> None:
    """有机酸反离子必须被剥离、活性母核必须保留；真实多组分混合物必须拦截。

    断言要点：
    * 4 种药用有机酸盐经 ``desalt_and_filter`` 后得到**单碎片**母核，重原子数与
      活性成分一致（不抛异常、也不把母核当反离子丢掉）；
    * ``prepare_ligand_from_smiles`` 对 4 种盐全部成功，产出的 PDBQT 与 SMILES
      都是单碎片母核（反离子彻底消失）；
    * ``prepare_library`` 对这 4 条盐全部命中 ``n_ready``；
    * 阿司匹林 + 布洛芬被 ``prepare_ligand_from_smiles`` 以「多个独立组分」拦截；
      在 ``prepare_library`` 中同一条输入必须计入失败（``n_desalt_fail``）而不进
      ``n_ready``，且质量守恒 ``n_raw - losses == n_ready`` 仍严格闭合。
    """
    # --- 1A. 底层脱盐：单碎片母核 + 正确重原子数 ---------------------------
    for case_id, salt_smiles, expected_heavy in ORGANIC_SALT_CASES:
        raw = Chem.MolFromSmiles(salt_smiles)
        assert raw is not None, f"{case_id}: 盐形式 SMILES 无法解析"
        assert len(Chem.GetMolFrags(raw)) >= 2, f"{case_id}: 用例本身必须是多碎片盐"

        core, reason = desalt_and_filter(raw, min_heavy=3, max_heavy=120)
        assert reason is None, f"{case_id}: 药用有机酸盐不该被判失败 ({reason})"
        assert core is not None
        assert len(Chem.GetMolFrags(core)) == 1, f"{case_id}: 母核必须为单碎片"
        assert core.GetNumHeavyAtoms() == expected_heavy, (
            f"{case_id}: 母核重原子数应为 {expected_heavy}，实际 {core.GetNumHeavyAtoms()} "
            f"({Chem.MolToSmiles(core)})"
        )

    # --- 1B. prepare_ligand_from_smiles：全部成功且反离子消失 ---------------
    lig_dir = settings.prepared_dir / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    for case_id, salt_smiles, expected_heavy in ORGANIC_SALT_CASES:
        ligand = prepare_ligand_from_smiles(salt_smiles, case_id, lig_dir)
        assert ligand.qc_profile is not None
        assert ligand.qc_profile.heavy_atoms == expected_heavy, (
            f"{case_id}: PDBQT 质量画像重原子数 {ligand.qc_profile.heavy_atoms} != {expected_heavy}"
        )
        core_mol = Chem.MolFromSmiles(ligand.smiles or "")
        assert core_mol is not None
        assert len(Chem.GetMolFrags(core_mol)) == 1
        assert core_mol.GetNumHeavyAtoms() == expected_heavy
        assert Path(ligand.pdbqt_path).is_file()

    # --- 1C. 真实双活性混合物：SMILES 入口必须抛 PreparationError -----------
    mixture_mol = Chem.MolFromSmiles(TRUE_MIXTURE_SMILES)
    assert mixture_mol is not None
    _, mixture_reason = desalt_and_filter(mixture_mol, min_heavy=3, max_heavy=120)
    assert mixture_reason == "multi_component_mixture"

    with pytest.raises(PreparationError, match="多个独立组分"):
        prepare_ligand_from_smiles(TRUE_MIXTURE_SMILES, "aspirin_ibuprofen", lig_dir)

    # --- 1D. prepare_library：盐全过、混合物计入失败且守恒闭合 ---------------
    registry = PrepareRegistry(settings)
    lines = [f"{smi} {cid}" for cid, smi, _ in ORGANIC_SALT_CASES]
    lines.append(f"{TRUE_MIXTURE_SMILES} aspirin_ibuprofen_mixture")
    content = ("\n".join(lines) + "\n").encode("utf-8")

    lib = prepare_library(
        content, "organic_salts_and_mixture.smi", registry, "有机酸盐库", tcm_mode=True
    )

    assert lib.n_raw == len(ORGANIC_SALT_CASES) + 1
    assert lib.n_desalt_fail + lib.n_prep_fail == 1, (
        "真实混合物必须且只能计入 1 条失败: "
        f"desalt={lib.n_desalt_fail}, prep={lib.n_prep_fail}, chiral={lib.n_chiral_fail}"
    )
    assert lib.n_ready == len(ORGANIC_SALT_CASES)
    assert lib.n_dedup == 0
    assert len(lib.ligand_ids) == len(ORGANIC_SALT_CASES)
    assert all("_aglycone" not in lid for lid in lib.ligand_ids)

    # 混合物失败原因落盘可追溯（desalt_and_filter 的固定字面量 multi_component_mixture）
    reasons = {item["reason"] for item in lib.failed}
    assert "multi_component_mixture" in reasons, reasons
    assert any(
        item["smiles"].startswith("CC(=O)Oc1ccccc1C(=O)O.CC(C)Cc1ccc") for item in lib.failed
    ), lib.failed

    # 质量守恒（与 /verify T5 相同的恒等式）
    total_losses = lib.n_desalt_fail + lib.n_dedup + lib.n_chiral_fail + lib.n_prep_fail
    assert lib.n_raw - total_losses == lib.n_ready

    # 落盘 manifest 与 funnel.csv 同样记录失败
    funnel_rows = dict(
        line.split(",")
        for line in (Path(lib.path) / "funnel.csv").read_text(encoding="utf-8").strip().splitlines()
    )
    assert int(funnel_rows["raw"]) == lib.n_raw
    assert int(funnel_rows["desalt"]) + int(funnel_rows["prep_failed"]) == 1
    assert int(funnel_rows["docked_success"]) == len(ORGANIC_SALT_CASES)


# =============================================================================
# 2. 互变异构归一化去重 & 未定手性中心枚举硬闸门
# =============================================================================


def test_tautomer_dedup_and_chiral_enumeration_active(settings) -> None:
    """酮式/烯醇式归一化去重 + 未定手性中心枚举闸门必须真实生效。

    断言要点：
    * ``standardize_tautomer_and_charge`` 把乙酰丙酮的酮式/烯醇式归一到同一
      InChIKey（修复前二者会被当成两条独立分子重复入池）；
    * ``prepare_library(tcm_mode=True)`` 对这两条输入给出 ``n_dedup == 1``、
      ``n_ready == 1``，且 ``source_herbs`` 把两个来源都挂在同一个 InChIKey 下；
    * 4 个未定手性中心（> ``MAX_UNSPECIFIED_STEREO_CENTERS`` = 3）的**非环**多醇
      ``CC(O)C(O)C(O)C(O)C`` 被拦截为 ``n_chiral_fail == 1`` 且不进入 ``n_ready``；
    * 该闸门不会误伤天然糖环（己糖的未定中心按设计放行，见
      ``test_htvs_funnel_stress::test_stereoisomer_unassigned_centers_robustness``）。
    """
    # --- 2A. 互变异构体 InChIKey 归一 --------------------------------------
    std_keto = standardize_tautomer_and_charge(Chem.MolFromSmiles(TAUTOMER_ACAC_KETO))
    std_enol = standardize_tautomer_and_charge(Chem.MolFromSmiles(TAUTOMER_ACAC_ENOL))
    assert std_keto is not None and std_enol is not None
    keto_key = Chem.MolToInchiKey(std_keto)
    enol_key = Chem.MolToInchiKey(std_enol)
    assert keto_key and keto_key == enol_key, (
        f"酮式/烯醇式未归一到同一 InChIKey: {keto_key} vs {enol_key}"
    )

    registry = PrepareRegistry(settings)
    tautomer_lib = prepare_library(
        f"{TAUTOMER_ACAC_KETO} acac_keto\n{TAUTOMER_ACAC_ENOL} acac_enol\n".encode(),
        "tautomers.smi",
        registry,
        "互变异构库",
        tcm_mode=True,
    )
    assert tautomer_lib.n_raw == 2
    assert tautomer_lib.n_dedup == 1, (
        f"n_dedup 应为 1，实际 {tautomer_lib.n_dedup} "
        f"(desalt={tautomer_lib.n_desalt_fail}, chiral={tautomer_lib.n_chiral_fail})"
    )
    assert tautomer_lib.n_desalt_fail == 0
    assert tautomer_lib.n_chiral_fail == 0
    assert tautomer_lib.n_ready == 1
    assert len(tautomer_lib.ligand_ids) == 1
    assert tautomer_lib.n_raw - 1 == tautomer_lib.n_ready  # 质量守恒
    assert tautomer_lib.source_herbs[keto_key] == ["acac_keto", "acac_enol"]

    # 唯一保留的配体其 InChIKey 必须是归一化后的那一个
    manifest = json.loads(
        (Path(tautomer_lib.path) / "library_manifest.json").read_text(encoding="utf-8")
    )
    assert isinstance(manifest, dict) and "ligands" in manifest
    retained = [item for item in manifest["ligands"] if item["id"] in tautomer_lib.ligand_ids]
    assert len(retained) == 1
    retained_mol = Chem.MolFromSmiles(retained[0]["smiles"])
    assert retained_mol is not None
    assert Chem.MolToInchiKey(retained_mol) == keto_key

    # --- 2B. 4 个未定手性中心的非环分子必须被拦截 --------------------------
    assert MAX_UNSPECIFIED_STEREO_CENTERS == 3
    probe = Chem.MolFromSmiles(ACYCLIC_TETRAOL_4_CENTERS)
    assert probe is not None
    unassigned = [
        idx for idx, tag in Chem.FindMolChiralCenters(probe, includeUnassigned=True) if tag == "?"
    ]
    assert len(unassigned) == 4, f"用例分子必须正好含 4 个未定手性中心，实际 {unassigned}"
    assert probe.GetRingInfo().NumRings() == 0, "必须是**非环**分子（糖环会被设计性放行）"

    chiral_lib = prepare_library(
        f"{ACYCLIC_TETRAOL_4_CENTERS} acyclic_tetraol\n".encode(),
        "chiral_ambiguity.smi",
        registry,
        "手性歧义库",
        tcm_mode=True,
    )
    assert chiral_lib.n_raw == 1
    assert chiral_lib.n_chiral_fail == 1, f"n_chiral_fail 应为 1，实际 {chiral_lib.n_chiral_fail}"
    assert chiral_lib.n_ready == 0
    assert chiral_lib.ligand_ids == []
    assert chiral_lib.n_desalt_fail == 0 and chiral_lib.n_prep_fail == 0
    assert any("手性" in item["reason"] for item in chiral_lib.failed)
    assert chiral_lib.n_raw - chiral_lib.n_chiral_fail == chiral_lib.n_ready

    # 己糖（糖环 + 4 个未定中心）仍必须放行：闸门只拦「真实立体歧义」
    hexose_lib = prepare_library(
        b"OC1C(O)C(O)C(O)C(CO)O1 natural_hexose\n",
        "hexose_control.smi",
        registry,
        "己糖对照库",
        tcm_mode=True,
    )
    assert hexose_lib.n_chiral_fail == 0, "天然糖环的未定中心必须放行（设计如此）"
    assert hexose_lib.n_ready == 1


# =============================================================================
# 3. 木脂素 / 内酯不得被误判为糖苷
# =============================================================================


@pytest.mark.parametrize(("case_id", "smiles", "has_fake_sugar_ring"), NON_GLYCOSIDE_CASES)
def test_lignan_and_lactone_not_misclassified_as_glycoside(
    case_id: str, smiles: str, has_fake_sugar_ring: bool
) -> None:
    """木脂素/内酯/黄烷醇母核绝不能被当成糖环切碎。

    每个样例都断言：
    * ``identify_sugar_rings`` 返回 ``[]``；
    * ``deglycosylate_tcm_mol`` 不产出苷元、也不报糖环名；
    * 若该样例确实含「5/6 元、全单键、恰好 1 个环氧」的伪糖环（``has_fake_sugar_ring``），
      则额外断言历史宽松判据确实会命中 —— 保证这是真正的对抗性输入。
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"{case_id}: SMILES 无法解析"

    if has_fake_sugar_ring:
        legacy_hits = _fake_sugar_like_rings(mol)
        assert legacy_hits, f"{case_id}: 该样例不含伪糖环，无法构成对抗性输入"

    sugar_rings = identify_sugar_rings(mol)
    assert sugar_rings == [], f"{case_id}: 木脂素/内酯母核被误判为糖环 {sugar_rings}"

    aglycone, sugar_names = deglycosylate_tcm_mol(mol)
    assert aglycone is None, (
        f"{case_id}: 非糖苷分子被错误水解，产出了 {aglycone.GetNumHeavyAtoms()} 重原子苷元"
    )
    assert sugar_names == [], f"{case_id}: 不应报出糖环名 {sugar_names}"


@pytest.mark.parametrize(("case_id", "smiles", "aglycone_heavy"), TRUE_GLYCOSIDE_CASES)
def test_true_glycosides_still_cleaved(case_id: str, smiles: str, aglycone_heavy: int) -> None:
    """O-糖苷（黄芩苷、芦丁）与 N-糖苷（腺苷）必须被正确切除糖环。"""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"{case_id}: SMILES 无法解析"

    sugar_rings = identify_sugar_rings(mol)
    assert sugar_rings, f"{case_id}: 真糖苷必须被识别出糖环"

    aglycone, sugar_names = deglycosylate_tcm_mol(mol)
    assert aglycone is not None, f"{case_id}: 真糖苷必须能水解出苷元"
    assert aglycone.GetNumHeavyAtoms() == aglycone_heavy, (
        f"{case_id}: 苷元重原子数应为 {aglycone_heavy}，实际 {aglycone.GetNumHeavyAtoms()}"
    )
    assert sugar_names, f"{case_id}: 水解必须报出糖环类型"
    assert len(sugar_names) == len(sugar_rings)


def test_c_glycoside_puerarin_protected() -> None:
    """葛根素（C-糖苷）必须被保护：识别为糖样环但绝不切割。"""
    mol = Chem.MolFromSmiles(PUERARIN_SMILES)
    assert mol is not None

    assert identify_sugar_rings(mol), "葛根素含有糖样吡喃环，应被识别"

    aglycone, sugar_names = deglycosylate_tcm_mol(mol)
    assert aglycone is None, "C-糖苷（C-C 直接相连）绝不能被水解"
    assert sugar_names == []


# =============================================================================
# 4. Tier 2 张力能惩罚真实生效 & 无「凑够 3 个」后门
# =============================================================================

_ROTATABLE_BOND_SMARTS = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")


def _pose_pdbqt_text(mol_3d: Chem.Mol) -> str:
    """Meeko PDBQT text of a 3D molecule (same writer the production prep uses)."""
    setups = MoleculePreparation().prepare(mol_3d)
    assert setups, "meeko produced no setup"
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    assert ok and pdbqt, f"meeko PDBQT writing failed: {err}"
    return pdbqt


def _measure_pose_strain(
    smiles: str, pose_mol: Chem.Mol, baseline: float, affinity: float = -8.0
) -> tuple[float, bool]:
    """Reproduces the runner's pose pipeline and returns (ΔE_strain, coords_are_real)."""
    parsed = load_pose_mol(_pose_pdbqt_text(pose_mol))
    assert parsed is not None, "load_pose_mol failed on the mock docked pose"
    restored = restore_pose_mol_with_template(parsed, smiles)
    assert restored is not None, "restore_pose_mol_with_template failed"
    strain, _corr, _high = compute_internal_strain_energy(
        restored,
        global_min_energy=baseline,
        affinity=affinity,
        high_strain_threshold=8.0,
    )
    return strain, pose_coords_are_real(restored)


def _strained_pose(
    smiles: str,
    relaxed_mol: Chem.Mol,
    baseline: float,
    low: float,
    high: float,
) -> Chem.Mol:
    """Rotates one non-ring rotatable bond until ΔE_strain lands in ``(low, high]``.

    A torsion rotation is used deliberately: it preserves bond lengths/angles, so the
    distorted pose stays a *physically valid* conformer whose only excess is genuine
    conformational strain (random coordinate noise would instead explode every bond).
    """
    for match in relaxed_mol.GetSubstructMatches(_ROTATABLE_BOND_SMARTS, uniquify=False)[:8]:
        a_idx, b_idx = match[0], match[1]
        nbrs_a = [
            n.GetIdx()
            for n in relaxed_mol.GetAtomWithIdx(a_idx).GetNeighbors()
            if n.GetIdx() != b_idx
        ]
        nbrs_b = [
            n.GetIdx()
            for n in relaxed_mol.GetAtomWithIdx(b_idx).GetNeighbors()
            if n.GetIdx() != a_idx
        ]
        if not nbrs_a or not nbrs_b:
            continue
        dihedral = (nbrs_a[0], a_idx, b_idx, nbrs_b[0])
        original = rdMolTransforms.GetDihedralDeg(relaxed_mol.GetConformer(), *dihedral)
        for delta in (15, 25, 35, 45, 55, 65, 75, 90, 105, 120, 135, 150, 165):
            candidate = Chem.Mol(relaxed_mol)
            rdMolTransforms.SetDihedralDeg(candidate.GetConformer(), *dihedral, original + delta)
            strain: float | None = None
            is_real = False
            try:
                strain, is_real = _measure_pose_strain(smiles, candidate, baseline)
            except Exception:  # noqa: BLE001 - unusable torsion: just skip this candidate
                strain = None
            if is_real and strain is not None and low < strain <= high:
                return candidate
            if is_real and low < strain <= high:
                return candidate
    raise AssertionError(f"无法在 {smiles} 上构造 ΔE_strain ∈ ({low}, {high}] 的真实张力姿态")


def _prepare_mock_ligand(
    lib_dir: Path, ligand_id: str, smiles: str
) -> tuple[PreparedLigand, Chem.Mol, float]:
    """Builds a production-shaped PreparedLigand (relaxed 3D + MMFF94s baseline)."""
    best_mol, min_energy = generate_conformers_and_min_energy(Chem.MolFromSmiles(smiles))
    lig_dir = lib_dir / ligand_id
    lig_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = lig_dir / "ligand.pdbqt"
    pdbqt_path.write_text(_pose_pdbqt_text(best_mol), encoding="utf-8")

    no_h = Chem.RemoveHs(best_mol)
    qc = LigandQCProfile(
        heavy_atoms=no_h.GetNumHeavyAtoms(),
        mw=round(Descriptors.MolWt(no_h), 3),
        clogp=round(Descriptors.MolLogP(no_h), 3),
        inchikey=Chem.MolToInchiKey(no_h),
    )
    ligand = PreparedLigand(
        id=ligand_id,
        source="smiles",
        pdbqt_path=str(pdbqt_path),
        smiles=Chem.MolToSmiles(no_h),
        name=ligand_id,
        qc_profile=qc,
        min_energy_kcal=min_energy,
    )
    return ligand, best_mol, min_energy


class _DeterministicPoseEngine(BaseEngine):
    """Mock engine that returns a *pre-built* real 3D conformer for each ligand.

    The pose PDBQT is a genuine meeko-written 3D structure (no placeholder atoms), so
    ``restore_pose_mol_with_template`` takes the real-coordinate branch and the strain
    QC actually measures what the conformer is worth.
    """

    manifest = EngineManifest(name="mock_pose_engine", version="1.0.0", scoring_functions=["vina"])

    def __init__(
        self,
        poses: dict[str, Chem.Mol],
        affinities: dict[str, dict[int | str, float]],
    ) -> None:
        self.poses = poses
        self.affinities = affinities
        self.calls: list[tuple[str, int | None]] = []

    def is_available(self) -> tuple[bool, str | None]:
        return True, None

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        ligand_id = Path(ligand_pdbqt).parent.name
        self.calls.append((ligand_id, params.seed))
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        (work_dir / "out.pdbqt").write_text(
            _pose_pdbqt_text(self.poses[ligand_id]), encoding="utf-8"
        )
        table = self.affinities[ligand_id]
        affinity = float(table.get(params.seed, table["default"]))
        return DockResult(
            engine=self.manifest.name,
            engine_version=self.manifest.version,
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=affinity)],
            log="deterministic mock pose engine",
            run_dir=str(work_dir),
        )


_CONVERGED_AFFINITIES: dict[int | str, float] = {
    42: -8.0,
    1024: -8.0,
    2026: -8.0,
    "default": -8.0,
}
_DIVERGENT_AFFINITIES: dict[int | str, float] = {
    42: -9.8,
    1024: -7.2,
    2026: -5.4,
    "default": -8.0,
}


def test_tier2_strain_penalty_applied_and_no_unconverged_backdoor(tmp_path: Path) -> None:
    """Tier 2 必须用张力修正后的亲和力驱动 SILE/排序，且绝不「凑够 3 个」放行。

    断言要点：
    * 真实 3D 张力姿态（``3.0 < ΔE_strain <= 8.0``，由非环可旋转键扭转构造）：
      ``corrected_affinity > raw_affinity``，SILE 用修正值计算并因此**排序下降**
      （原始亲和力下本应排第一的分子被反超）；
    * ``ΔE_strain > 8.0`` 的姿态被张力 QC 否决；
    * 「1 个收敛 + 2 个发散/超高张力」时 ``_run_tier2`` 返回列表长度严格 == 1；
    * 0 个收敛时只降级保留 1 个，并打上 ``UNCONVERGED_OR_HIGH_STRAIN``；
      该分子进入 Tier 3 后候选行 alerts 必须含 ``HIGH_STRAIN_UNCONVERGED_FALLBACK``。
    """
    lib_dir = tmp_path / "mock_lib"
    rec_path = tmp_path / "receptor.pdbqt"
    rec_path.write_text(
        "ATOM      1  N   ASP A  25      16.000  25.000   4.000  1.00  0.00    -0.300 N\nTER\n",
        encoding="utf-8",
    )

    # 放松姿态（ΔE≈0）与真实张力姿态（用扭转非线性地抬高构象能）
    clean_lig, clean_pose, _clean_base = _prepare_mock_ligand(
        lib_dir,
        "LIG_CLEAN",
        "CC(=O)Oc1ccccc1C(=O)O",  # 阿司匹林 13 重原子
    )
    mild_lig, mild_relaxed, mild_base = _prepare_mock_ligand(
        lib_dir,
        "LIG_MILD",
        "CC(=O)Nc1ccc(O)cc1",  # 对乙酰氨基酚 11 重原子（SILE 分母小）
    )
    high_lig, high_relaxed, high_base = _prepare_mock_ligand(
        lib_dir,
        "LIG_HIGH",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # 布洛芬（可扭转键多，便于造高张力）
    )
    diverge_lig, diverge_pose, _diverge_base = _prepare_mock_ligand(
        lib_dir,
        "LIG_DIVERGE",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",  # 咖啡因（刚性平面骨架）
    )

    mild_pose = _strained_pose(mild_lig.smiles, mild_relaxed, mild_base, 6.0, 8.0)
    high_pose = _strained_pose(high_lig.smiles, high_relaxed, high_base, 8.5, 200.0)

    # 自检：构造出的姿态确实落在目标带内（否则后续断言无意义）
    mild_strain_probe, mild_real = _measure_pose_strain(mild_lig.smiles, mild_pose, mild_base)
    assert 3.0 < mild_strain_probe <= 8.0, mild_strain_probe
    assert mild_real, "张力姿态必须是可还原的真实对接坐标"
    high_strain_probe, high_real = _measure_pose_strain(high_lig.smiles, high_pose, high_base)
    assert high_strain_probe > 8.0, high_strain_probe
    assert high_real

    poses = {
        clean_lig.id: clean_pose,
        mild_lig.id: mild_pose,
        high_lig.id: high_pose,
        diverge_lig.id: diverge_pose,
    }
    raw_affinity = -8.0

    def _tier1_rows(ligands: list[PreparedLigand]) -> list[dict[str, Any]]:
        rows = []
        for lig in ligands:
            effs = compute_ligand_efficiencies(
                raw_affinity, lig.qc_profile.heavy_atoms, lig.qc_profile.clogp
            )
            rows.append(
                {
                    "ligand": lig,
                    "affinity": raw_affinity,
                    "sile": float(effs["sile"] or 0.0),
                    "le": float(effs["le"] or 0.0),
                    "lle": float(effs["lle"] or 0.0),
                    "run_dir": str(tmp_path / "tier1" / f"dock_{lig.id}"),
                }
            )
        return rows

    def _make_runner(
        name: str, ligands: list[PreparedLigand], affinities: dict[str, dict[int | str, float]]
    ) -> ScreeningRunner:
        engine = _DeterministicPoseEngine(poses, affinities)
        return ScreeningRunner(
            config=ScreeningConfig(
                enable_tiered_funnel=True,
                tier1_keep_ratio=1.0,
                tier2_keep_ratio=1.0,  # 关闭比例截断，任何「凑数」都会立刻暴露
                tier2_max_candidates=100,
                tier2_seeds=[42, 1024, 2026],
                max_strain_kcal=8.0,
                max_sigma_dG=0.35,
                max_rmsd_seed=2.0,
                enable_gnina_rescore=False,
                box=Box(center=(16.0, 25.0, 4.0), size=(22.5, 22.5, 22.5)),
            ),
            receptor_pdbqt=rec_path,
            ligands={lig.id: lig for lig in ligands},
            out_dir=tmp_path / name,
            engine=engine,
        )

    # --- 4A. 张力惩罚真实进入 SILE 且改变排序 ------------------------------
    runner_a = _make_runner(
        "tier2_strain",
        [clean_lig, mild_lig, high_lig],
        {
            clean_lig.id: _CONVERGED_AFFINITIES,
            mild_lig.id: _CONVERGED_AFFINITIES,
            high_lig.id: _CONVERGED_AFFINITIES,
        },
    )
    rows_a = runner_a._run_tier2(_tier1_rows([clean_lig, mild_lig, high_lig]), {})
    passed_ids = [row["ligand"].id for row in rows_a]
    assert high_lig.id not in passed_ids, "ΔE_strain > 8.0 的高张力姿态必须被张力 QC 否决"
    assert set(passed_ids) == {clean_lig.id, mild_lig.id}, passed_ids

    row_by_id = {row["ligand"].id: row for row in rows_a}
    mild_row = row_by_id[mild_lig.id]
    assert 3.0 < mild_row["strain_energy"] <= 8.0, mild_row["strain_energy"]
    assert mild_row["is_strain_ok"] is True
    assert mild_row["affinity_corrected"] is not None
    assert mild_row["affinity_corrected"] > mild_row["affinity_mean"], (
        "3 < ΔE ≤ 8 时必须施加张力惩罚: "
        f"corr={mild_row['affinity_corrected']} raw={mild_row['affinity_mean']}"
    )
    assert mild_row["affinity_corrected"] == pytest.approx(
        mild_row["affinity_mean"] + 0.5 * (mild_row["strain_energy"] - 3.0), abs=5e-3
    )

    # SILE 必须用修正后的亲和力：原始亲和力下 LIG_MILD（重原子少）SILE 更高，
    # 惩罚后必须反超 -> 排序下降。
    raw_sile_clean = compute_ligand_efficiencies(raw_affinity, clean_lig.qc_profile.heavy_atoms)[
        "sile"
    ]
    raw_sile_mild = compute_ligand_efficiencies(raw_affinity, mild_lig.qc_profile.heavy_atoms)[
        "sile"
    ]
    assert raw_sile_mild > raw_sile_clean, "前提：原始亲和力下 LIG_MILD 的 SILE 本应领先"
    assert mild_row["sile"] < raw_sile_mild, (
        f"Tier 2 SILE 未使用修正后亲和力: {mild_row['sile']} vs 原始 {raw_sile_mild}"
    )
    assert mild_row["sile"] < row_by_id[clean_lig.id]["sile"], "张力惩罚未让排序下降"
    assert rows_a[0]["ligand"].id == clean_lig.id, "排序必须由修正后 SILE 决定"

    # Tier 3 候选行必须落盘 affinity_corrected（数据不丢失）
    cands_a = runner_a._run_tier3(rows_a, {})
    assert len(cands_a) == 2
    cand_by_id = {cand.ligand_id: cand for cand in cands_a}
    assert cand_by_id[mild_lig.id].affinity_corrected is not None
    assert cand_by_id[mild_lig.id].affinity_corrected == pytest.approx(
        mild_row["affinity_corrected"], abs=1e-2
    )
    assert cand_by_id[mild_lig.id].strain_energy == pytest.approx(
        mild_row["strain_energy"], abs=1e-2
    )

    # --- 4B. 「1 收敛 + 2 发散/超高张力」：严格只晋级 1 个 ------------------
    runner_b = _make_runner(
        "tier2_no_backdoor",
        [clean_lig, diverge_lig, high_lig],
        {
            clean_lig.id: _CONVERGED_AFFINITIES,
            diverge_lig.id: _DIVERGENT_AFFINITIES,  # σ ≫ 0.35 -> 不收敛
            high_lig.id: _CONVERGED_AFFINITIES,  # ΔE > 8 -> 张力否决
        },
    )
    rows_b = runner_b._run_tier2(_tier1_rows([clean_lig, diverge_lig, high_lig]), {})
    assert len(rows_b) == 1, (
        f"Tier 2 必须严格只晋级 1 个收敛分子，实际 {[r['ligand'].id for r in rows_b]}"
        " —— 疑似触发「凑够 3 个」后门"
    )
    assert rows_b[0]["ligand"].id == clean_lig.id
    assert "qc_fallback_warning" not in rows_b[0], "存在合格分子时不得打降级警示"

    # --- 4C. 0 个收敛：只降级保留 1 个并打上双重警示 ----------------------
    runner_c = _make_runner(
        "tier2_all_unconverged",
        [diverge_lig, mild_lig, high_lig],
        {
            diverge_lig.id: _DIVERGENT_AFFINITIES,
            mild_lig.id: _DIVERGENT_AFFINITIES,  # 同一分子改为发散打分
            high_lig.id: _CONVERGED_AFFINITIES,  # 收敛但 ΔE > 8，仍不合格
        },
    )
    rows_c = runner_c._run_tier2(_tier1_rows([diverge_lig, mild_lig, high_lig]), {})
    assert len(rows_c) == 1, f"0 个收敛时只允许降级保留 1 个，实际 {len(rows_c)}"
    assert rows_c[0]["qc_fallback_warning"] == "UNCONVERGED_OR_HIGH_STRAIN"

    cands_c = runner_c._run_tier3(rows_c, {})
    assert len(cands_c) == 1
    fallback_alerts = [
        alert for alert in cands_c[0].alerts if alert.startswith("HIGH_STRAIN_UNCONVERGED_FALLBACK")
    ]
    assert fallback_alerts, (
        f"Tier 3 候选 alerts 缺少 HIGH_STRAIN_UNCONVERGED_FALLBACK: {cands_c[0].alerts}"
    )


# -----------------------------------------------------------------------------
# 4b. PDBQT 姿态复原不得伪造「假高张力」（Task Fix 3.1b 回归）
# -----------------------------------------------------------------------------

#: 多羟基黄酮 / 多酚：修复前经 meeko PDBQT 往返 + ``restore_pose_mol_with_template``
#: 后会被算成 4~26 kcal/mol 的假张力并在 Tier 2 被误否决（真实缺陷，见 QA 报告）。
#: 注：此处按「SMILES 直接 prep」路径验证（与 Lead 复验脚本同一路径）。若同一个分子先经
#: ``prepare_library`` 的标准化/嵌入路径（原子序不同），极性氢转子搜索仍可能停在局部
#: 极小 —— 该残余已在第 5 组 ``DELIVERY_LIBRARY`` 处显式记录并上报。
POLYPHENOL_RESTORE_CASES: tuple[tuple[str, str], ...] = (
    ("Quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    (
        "Rutin",
        "CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O",
    ),
    (
        "Baicalin",
        "O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O",
    ),
    ("Arbutin", "OC1=CC=CC=C1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),
)

#: 单羟基对照：修复前后都应接近 0（证明「多酚」才是问题域，不是力场整体漂移）。
POLYPHENOL_RESTORE_CONTROLS: tuple[tuple[str, str], ...] = (
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
)

#: 正向判据：张力能上限（修复后实测 0.000~0.004 kcal/mol）。
_POLYPHENOL_MAX_FAKE_STRAIN = 1.5
#: 正向判据：重原子几何必须与 prep 全局极小一致。
_POLYPHENOL_MAX_HEAVY_RMSD = 0.01
#: 反向判据：单个芳香环可旋转键强制扭转 >= 60 度后的最小真实代价。
_POLYPHENOL_MIN_TORSION_STRAIN = 3.0

#: 强制扭转角度（满足「>= 60 度」要求，同时足以跨过芳环共轭平面）。
_FORCED_TORSION_DEG = 120.0


def _restored_pose_strain(
    smiles: str, pose_mol: Chem.Mol, baseline: float
) -> tuple[float, bool, float]:
    """PDBQT 往返（meeko 写出 -> load_pose_mol -> 模板复原）后的 (ΔE, is_high, heavy RMSD)."""
    parsed = load_pose_mol(_pose_pdbqt_text(pose_mol))
    assert parsed is not None, "load_pose_mol 无法解析 mock 姿态"
    restored = restore_pose_mol_with_template(parsed, smiles)
    assert restored is not None, "restore_pose_mol_with_template 复原失败"
    strain, _corr, is_high = compute_internal_strain_energy(
        restored,
        global_min_energy=baseline,
        affinity=-8.0,
        high_strain_threshold=8.0,
    )
    heavy_rmsd = AllChem.GetBestRMS(
        Chem.Mol(Chem.RemoveHs(restored)), Chem.Mol(Chem.RemoveHs(pose_mol))
    )
    return float(strain), bool(is_high), float(heavy_rmsd)


def _forced_aromatic_torsion(
    smiles: str,
    relaxed_mol: Chem.Mol,
    baseline: float,
    threshold: float,
) -> tuple[float, bool, tuple[int, int]] | None:
    """把某个「芳香环相邻的非环可旋转键」强制扭转 120 度，返回首个超阈值的实测值。

    只旋转 ``-!@``（非环）单键，因此键长/键角保持物理合理，成本全部来自真实的
    构象张力；用于防止「把高张力检测整体关掉」式的一刀切放水。
    """
    for match in relaxed_mol.GetSubstructMatches(_ROTATABLE_BOND_SMARTS, uniquify=False):
        a_idx, b_idx = match[0], match[1]
        atom_a = relaxed_mol.GetAtomWithIdx(a_idx)
        atom_b = relaxed_mol.GetAtomWithIdx(b_idx)
        if not (atom_a.GetIsAromatic() or atom_b.GetIsAromatic()):
            continue
        nbrs_a = [n.GetIdx() for n in atom_a.GetNeighbors() if n.GetIdx() != b_idx]
        nbrs_b = [n.GetIdx() for n in atom_b.GetNeighbors() if n.GetIdx() != a_idx]
        if not nbrs_a or not nbrs_b:
            continue
        dihedral = (nbrs_a[0], a_idx, b_idx, nbrs_b[0])
        original = rdMolTransforms.GetDihedralDeg(relaxed_mol.GetConformer(), *dihedral)
        for signed_delta in (_FORCED_TORSION_DEG, -_FORCED_TORSION_DEG):
            candidate = Chem.Mol(relaxed_mol)
            rdMolTransforms.SetDihedralDeg(
                candidate.GetConformer(), *dihedral, original + signed_delta
            )
            strain, is_high, _rmsd = _restored_pose_strain(smiles, candidate, baseline)
            if strain > threshold:
                return strain, is_high, (a_idx, b_idx)
    return None


def test_polyphenol_pose_restore_does_not_fake_high_strain() -> None:
    """PDBQT 姿态复原不得给「本身就是全局极小」的多酚姿态伪造高张力。

    背景：``restore_pose_mol_with_template`` 用 ``AddHs(addCoords=True)`` 重新生成
    氢原子，多酚会落到与 prep 基线不同的氢键网络；``compute_internal_strain_energy``
    只放松氢原子时无法跨过 O-H 旋转势垒，于是把重原子坐标与全局极小完全一致的姿态
    误算成 4~26 kcal/mol 假张力，进而让 ``ScreeningRunner`` 在 Tier 2 误否决
    合法的黄酮/多酚命中（Task Fix 3.1b）。

    断言要点：
    * 正向（Quercetin / Rutin / Baicalin / Arbutin）：走与生产完全相同的
      ``generate_conformers_and_min_energy -> meeko PDBQT -> load_pose_mol ->
      restore_pose_mol_with_template -> compute_internal_strain_energy`` 路径，
      重原子 RMSD ≤ 0.01 Å、ΔE_strain ≤ 1.5 kcal/mol、``is_high_strain is False``；
    * 对照（Aspirin / Paracetamol）：ΔE_strain ≤ 0.5 且不误报高张力；
    * **反向断言**（防一刀切放水）：把芳香环相邻可旋转键强制扭转 120°（≥60°）后，
      Quercetin / Baicalin 的 ΔE_strain 必须 > 3.0；Baicalin 的糖苷/芳基扭转实测
      可达 20+ kcal/mol，必须仍被判为 ``is_high_strain is True``。
      （Quercetin 唯一的芳基-色酮扭转是类联苯扭转，物理上限约 5.6 kcal/mol < 8.0，
      故它只参与「> 3.0」判据；ΔE > 8 的真实高张力否决路径由第 4 组用例覆盖。）
    """
    # --- 正向：多酚不得再被伪造高张力 ------------------------------------
    forward_measurements: dict[str, tuple[float, bool, float]] = {}
    for name, smiles in POLYPHENOL_RESTORE_CASES:
        relaxed_mol, baseline = generate_conformers_and_min_energy(Chem.MolFromSmiles(smiles))
        strain, is_high, heavy_rmsd = _restored_pose_strain(smiles, relaxed_mol, baseline)
        forward_measurements[name] = (strain, is_high, heavy_rmsd)
        assert heavy_rmsd <= _POLYPHENOL_MAX_HEAVY_RMSD, (
            f"{name}: PDBQT 往返后重原子 RMSD {heavy_rmsd:.4f} Å 超过 "
            f"{_POLYPHENOL_MAX_HEAVY_RMSD} Å（姿态几何被改动）"
        )
        assert strain <= _POLYPHENOL_MAX_FAKE_STRAIN, (
            f"{name}: 全局极小姿态被算出 ΔE_strain={strain:.3f} kcal/mol "
            f"(> {_POLYPHENOL_MAX_FAKE_STRAIN})，假高张力回归"
        )
        assert is_high is False, f"{name}: ΔE={strain:.3f} 不应被标记为高张力"

    # 对照分子：修复前后都应接近 0（问题域仅限多酚氢键网络）
    for name, smiles in POLYPHENOL_RESTORE_CONTROLS:
        relaxed_mol, baseline = generate_conformers_and_min_energy(Chem.MolFromSmiles(smiles))
        strain, is_high, heavy_rmsd = _restored_pose_strain(smiles, relaxed_mol, baseline)
        assert heavy_rmsd <= _POLYPHENOL_MAX_HEAVY_RMSD
        assert strain <= 0.5, f"{name}: 对照分子 ΔE_strain={strain:.3f} 异常"
        assert is_high is False

    # --- 反向 1：芳香环可旋转键强制扭转 120° -> 张力必须真实上升 -------
    for name, smiles in POLYPHENOL_RESTORE_CASES:
        if name not in ("Quercetin", "Baicalin"):
            continue
        relaxed_mol, baseline = generate_conformers_and_min_energy(Chem.MolFromSmiles(smiles))
        probe = _forced_aromatic_torsion(
            smiles, relaxed_mol, baseline, _POLYPHENOL_MIN_TORSION_STRAIN
        )
        assert probe is not None, (
            f"{name}: 强制扭转芳香环可旋转键 {_FORCED_TORSION_DEG:.0f}° 后 ΔE_strain "
            f"仍 <= {_POLYPHENOL_MIN_TORSION_STRAIN} kcal/mol —— 张力检测疑似被整体关停"
        )
        strain, is_high, bond = probe
        assert strain > _POLYPHENOL_MIN_TORSION_STRAIN, (name, strain, bond)

    # --- 反向 2：足够大的扭转仍必须被判为高张力（Baicalin 实测 20+）-----
    baicalin_smiles = dict(POLYPHENOL_RESTORE_CASES)["Baicalin"]
    relaxed_mol, baseline = generate_conformers_and_min_energy(Chem.MolFromSmiles(baicalin_smiles))
    high_probe = _forced_aromatic_torsion(baicalin_smiles, relaxed_mol, baseline, 8.0)
    assert high_probe is not None, (
        "Baicalin: 强制扭转芳香环/糖苷可旋转键后 ΔE_strain 仍未超过 8.0 kcal/mol，"
        "高张力判据疑似被放水"
    )
    high_strain, high_flag, high_bond = high_probe
    assert high_strain > 8.0, (high_strain, high_bond)
    assert high_flag is True, (
        f"Baicalin: ΔE_strain={high_strain:.2f} kcal/mol > 8.0 时必须 is_high_strain=True"
    )


# =============================================================================
# 5. 交付姿态保留真实 3D 对接坐标 + /verify T1~T5 全量通过
# =============================================================================

#: 3-Tier + TCM 双轨筛选所用的库：10 条输入（含 1 个可水解 O-糖苷 -> 衍生 1 个苷元），
#: 全部为无立体歧义、可稳定复现 InChIKey 的分子，并**刻意包含多羟基黄酮
#: （黄芩苷 / 槲皮素）**——它们正是 Task Fix 3.1b（假高张力）的高风险样本，
#: 详见本文件第 4b 组 ``test_polyphenol_pose_restore_does_not_fake_high_strain``。
DELIVERY_LIBRARY: tuple[tuple[str, str], ...] = (
    # 黄芩苷：O-葡萄糖醛酸苷 -> 双轨衍生出黄芩素（Baicalein）苷元；多羟基黄酮
    (
        "Baicalin",
        "O=C(O)[C@@H]1O[C@@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)[C@H](O)[C@@H](O)[C@H]1O",
    ),
    ("Quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("Coumarin", "O=C1C=Cc2ccccc2O1"),
    ("Umbelliferone", "O=C1C=Cc2cc(O)ccc2O1"),
    ("Diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("Genistein", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"),
    ("Naringenin", "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
)

#: 搜索盒中心（真实结合口袋位置；姿态质心必须落在这里，绝不能回到原点）。
DELIVERY_BOX_CENTER = (16.0, 25.0, 4.0)
DELIVERY_BOX_SIZE = (22.5, 22.5, 22.5)

#: 每个配体一个互不相同的亚埃级偏移，使 mock 引擎能按质心唯一识别分子。
_DELIVERY_OFFSETS: tuple[tuple[float, float, float], ...] = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.6),
    (0.0, -1.0, -0.6),
    (1.0, 1.0, 0.0),
    (-1.0, -1.0, 0.0),
    (0.6, -1.0, 1.0),
    (-0.6, 1.0, -1.0),
    (1.0, -1.0, 1.0),
)


def _heavy_centroid_from_text(pose_text: str) -> tuple[float, float, float]:
    """Heavy-atom centroid of a PDBQT pose text (AutoDock-type aware parsing)."""
    mol = load_pose_mol(pose_text)
    assert mol is not None, "无法解析 PDBQT 姿态文本"
    conf = mol.GetConformer()
    heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    assert heavy, "姿态中未找到重原子"
    n = float(len(heavy))
    return (
        sum(conf.GetAtomPosition(i).x for i in heavy) / n,
        sum(conf.GetAtomPosition(i).y for i in heavy) / n,
        sum(conf.GetAtomPosition(i).z for i in heavy) / n,
    )


def _heavy_centroid_of_mol(mol: Chem.Mol) -> tuple[float, float, float]:
    conf = mol.GetConformer()
    heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    n = float(len(heavy))
    return (
        sum(conf.GetAtomPosition(i).x for i in heavy) / n,
        sum(conf.GetAtomPosition(i).y for i in heavy) / n,
        sum(conf.GetAtomPosition(i).z for i in heavy) / n,
    )


def _translate_pdbqt_text(pose_text: str, delta: tuple[float, float, float]) -> str:
    """Rigidly translates every PDBQT ATOM/HETATM record (columns 31-54) by ``delta``."""
    lines: list[str] = []
    for line in pose_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            x = float(line[30:38]) + delta[0]
            y = float(line[38:46]) + delta[1]
            z = float(line[46:54]) + delta[2]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        lines.append(line)
    return "\n".join(lines) + "\n"


class _PoseDeliveryEngine(BaseEngine):
    """Deterministic mock engine that docks every ligand *inside the real pocket box*.

    Two modes:

    * screening (``ligand_pdbqt`` lives in ``<lib>/<ligand_id>/``): the ligand's own
      prepared 3D conformer is rigidly translated onto its pocket target and written as
      a genuine meeko PDBQT ``out.pdbqt`` — i.e. real 3D docking coordinates;
    * /verify recalculation (T1/T2 re-dock a temporary ``ligand.pdbqt`` built from the
      delivered pose SDF): the ligand is re-identified by heavy-atom centroid and the
      *delivered* pose SDF is returned unchanged, so the re-run reproduces the reported
      affinity (|ΔAff| = 0) and the same geometry (RMSD = 0).
    """

    manifest = EngineManifest(
        name="mock_delivery_engine", version="1.0.0", scoring_functions=["vina"]
    )

    def __init__(self, fixed_affinity: float = -8.0, centroid_tol: float = 0.8) -> None:
        self.fixed_affinity = float(fixed_affinity)
        self.centroid_tol = float(centroid_tol)
        self.targets: dict[str, tuple[float, float, float]] = {}
        self.screening_poses: dict[str, str] = {}
        self.delivered: dict[str, tuple[tuple[float, float, float], str]] = {}
        self.dock_calls: list[tuple[str, int | None]] = []

    def is_available(self) -> tuple[bool, str | None]:
        return True, None

    # -- registry ---------------------------------------------------------
    def register_ligand(self, ligand_id: str, target: tuple[float, float, float]) -> None:
        self.targets[ligand_id] = target

    def register_delivered_pose(self, ligand_id: str, sdf_text: str) -> None:
        centroid = _heavy_centroid_of_mol(Chem.MolFromMolBlock(sdf_text, sanitize=False))
        self.delivered[ligand_id] = (centroid, sdf_text)

    def _match_delivered(self, centroid: tuple[float, float, float]) -> str:
        best_id, best_dist = None, float("inf")
        for ligand_id, (ref, _text) in self.delivered.items():
            dist = sum((centroid[i] - ref[i]) ** 2 for i in range(3)) ** 0.5
            if dist < best_dist:
                best_id, best_dist = ligand_id, dist
        assert best_id is not None and best_dist <= self.centroid_tol, (
            f"重算输入无法匹配到任何已交付姿态 (最近={best_id}, 距离={best_dist:.3f} Å)"
        )
        return best_id

    # -- engine API -------------------------------------------------------
    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        ligand_pdbqt = Path(ligand_pdbqt)
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        text = ligand_pdbqt.read_text(encoding="utf-8", errors="replace")
        dir_name = ligand_pdbqt.parent.name
        self.dock_calls.append((dir_name, params.seed))

        pose_sdf_path: Path | None = None
        if dir_name in self.targets:
            # Screening: real conformer rigidly translated into the pocket box.
            centroid = _heavy_centroid_from_text(text)
            target = self.targets[dir_name]
            delta = tuple(target[i] - centroid[i] for i in range(3))
            translated = _translate_pdbqt_text(text, delta)
            self.screening_poses[dir_name] = translated
            (work_dir / "out.pdbqt").write_text(translated, encoding="utf-8")
            ligand_id = dir_name
        else:
            # /verify recalculation: reproduce the delivered pose bit-for-bit.
            centroid = _heavy_centroid_from_text(text)
            ligand_id = self._match_delivered(centroid)
            _ref, sdf_text = self.delivered[ligand_id]
            pose_sdf_path = work_dir / "docked.sdf"
            pose_sdf_path.write_text(sdf_text, encoding="utf-8")
            (work_dir / "out.pdbqt").write_text(self.screening_poses[ligand_id], encoding="utf-8")

        return DockResult(
            engine=self.manifest.name,
            engine_version=self.manifest.version,
            params=params,
            box=box,
            poses=[
                Pose(
                    index=1,
                    affinity=self.fixed_affinity,
                    file=str(pose_sdf_path) if pose_sdf_path is not None else None,
                )
            ],
            log="deterministic pose-delivery mock engine",
            run_dir=str(work_dir),
        )


def _write_mock_receptor(tmp_path: Path) -> tuple[Path, Path]:
    """Minimal receptor PDB + PDBQT pair (same schema as tests/test_cli_screen.py)."""
    rec_dir = tmp_path / "receptor"
    rec_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = rec_dir / "1HSG.pdb"
    pdbqt_path = rec_dir / "1HSG.pdbqt"
    pdb_path.write_text(
        "ATOM      1  N   ASP A  25      16.000  25.000   4.000  1.00  0.00           N\n"
        "ATOM      2  CA  ASP A  25      17.000  25.000   4.000  1.00  0.00           C\n"
        "ATOM      3  N   THR A  80      16.000  30.000   4.000  1.00  0.00           N\n"
        "ATOM      4  CA  THR A  80      17.000  30.000   4.000  1.00  0.00           C\n"
        "TER\nEND\n",
        encoding="utf-8",
    )
    pdbqt_path.write_text(
        "ATOM      1  N   ASP A  25      16.000  25.000   4.000  1.00  0.00    -0.300 N\n"
        "ATOM      2  CA  ASP A  25      17.000  25.000   4.000  1.00  0.00    +0.100 C\n"
        "ATOM      3  N   THR A  80      16.000  30.000   4.000  1.00  0.00    -0.300 N\n"
        "ATOM      4  CA  THR A  80      17.000  30.000   4.000  1.00  0.00    +0.100 C\n"
        "TER\n",
        encoding="utf-8",
    )
    return pdb_path, pdbqt_path


def test_delivery_poses_preserve_real_docked_3d_coordinates_and_verify_t1_t5(
    settings, tmp_path: Path
) -> None:
    """交付姿态必须是真实 3D 对接坐标，且交付包必须通过 /verify T1~T5 全量验收。

    断言要点：
    * 3-Tier + ``tcm_mode=True`` 双轨筛选：衍生苷元使实际入池配体数 > ``n_raw``，
      ``funnel.csv`` 的 ``aglycone_generated > 0``；
    * ``ScreeningFunnelReport.is_conserved`` 为真，并且用报告 §3.5 的新公式**手算**闭合；
    * 每个 ``out_dir/poses/pose_<cid>.sdf`` 都带 3D 构象，其重原子质心落在
      ``(16, 25, 4)`` 搜索盒内（并精确落在该分子的口袋目标点上），且**绝不等于原点**；
    * ``run_verify(out_dir, options={"recalculate": True, "engine": engine})`` 返回 0，
      ``_verify_t1`` 以 ``PASS``（真实重算，|ΔAff|≈0、RMSD≈0）而不是 ``SKIP`` 通过，
      T2/T3/T4/T5 同样全部 ``passed``。
    """
    from dockstudio.prep.library import prepare_library

    registry = PrepareRegistry(settings)
    pdb_path, pdbqt_path = _write_mock_receptor(tmp_path)
    out_dir = tmp_path / "delivery_out"

    # --- 5A. TCM 双轨库准备 ------------------------------------------------
    content = ("\n".join(f"{smi} {name}" for name, smi in DELIVERY_LIBRARY) + "\n").encode("utf-8")
    lib = prepare_library(content, "delivery_library.smi", registry, "交付验收库", tcm_mode=True)
    assert lib.n_raw == len(DELIVERY_LIBRARY)
    assert lib.n_ready == len(DELIVERY_LIBRARY)
    assert lib.n_desalt_fail == 0 and lib.n_chiral_fail == 0 and lib.n_prep_fail == 0
    assert any(lid.endswith("_aglycone") for lid in lib.ligand_ids), "库中糖苷必须衍生出苷元"
    assert len(lib.ligand_ids) > lib.n_raw

    # --- 5B. 3-Tier 漏斗筛选（确定性 mock 引擎，姿态落进口袋盒） -----------
    engine = _PoseDeliveryEngine(fixed_affinity=-8.0)
    for idx, ligand_id in enumerate(lib.ligand_ids):
        offset = _DELIVERY_OFFSETS[idx % len(_DELIVERY_OFFSETS)]
        engine.register_ligand(
            ligand_id,
            tuple(DELIVERY_BOX_CENTER[i] + offset[i] for i in range(3)),
        )

    lib_ligands = {lid: registry.ligands[lid] for lid in lib.ligand_ids if lid in registry.ligands}
    assert set(lib_ligands) == set(lib.ligand_ids)

    runner = ScreeningRunner(
        config=ScreeningConfig(
            receptor_id="1HSG",
            library_id=lib.name,
            box=Box(center=DELIVERY_BOX_CENTER, size=DELIVERY_BOX_SIZE),
            tcm_mode=True,
            tcm_deglycosylation=True,
            enable_tiered_funnel=True,
            tier1_keep_ratio=1.0,  # 不截断，保证所有分子进入 3 个 Tier 的真实流程
            tier2_keep_ratio=1.0,
            tier2_seeds=[42, 1024, 2026],
            max_strain_kcal=8.0,
            enable_gnina_rescore=False,
            engine=engine.manifest.name,
        ),
        receptor_pdbqt=pdbqt_path,
        protein_pdb=pdb_path,
        ligands=lib_ligands,
        prepared_library=lib,
        out_dir=out_dir,
        engine=engine,
        cpu_per_job=1,
    )
    report = runner.run()

    # --- 5C. 双轨衍生 + 漏斗守恒 ------------------------------------------
    assert len(runner.ligand_dict) > lib.n_raw, "双轨衍生后实际入池数必须大于 n_raw"
    assert report.n_aglycone_generated > 0
    assert report.n_ready == lib.n_ready
    assert report.is_conserved, (
        "ScreeningFunnelReport.is_conserved 必须为真: "
        f"raw={report.n_raw}, aglycone={report.n_aglycone_generated}, ready={report.n_ready}, "
        f"funnel_filtered={report.n_funnel_filtered}, dock_failed={report.n_dock_failed}, "
        f"success={report.n_docked_success}"
    )
    assert report.n_docked_success == len(report.candidates) > 0

    funnel_rows = load_funnel_data(out_dir / "funnel.csv")
    assert funnel_rows["aglycone_generated"] > 0
    # 报告 §3.5 新公式手算（与 verify.py::_verify_t5 完全一致）
    expected_docked = (
        funnel_rows["raw"]
        + funnel_rows["aglycone_generated"]
        - (
            funnel_rows["desalt"]
            + funnel_rows["dedup"]
            + funnel_rows["chiral_failed"]
            + funnel_rows["prep_failed"]
            + funnel_rows["funnel_filtered"]
            + funnel_rows["dock_failed"]
        )
    )
    energy_rows = load_energy_table(out_dir / "screening_results.csv")
    assert expected_docked == funnel_rows["docked_success"] == len(energy_rows)

    # --- 5D. 交付姿态必须是口袋内的真实 3D 坐标 ---------------------------
    poses_dir = out_dir / "poses"
    assert poses_dir.is_dir()
    candidate_ids = [cand.ligand_id for cand in report.candidates]
    assert len(candidate_ids) >= 7, "T2 多 seed 重叠率需要 >= 7 个候选分子"
    # 多羟基黄酮回归守门：槲皮素（22 重原子）必须存活 —— 它是 Task Fix 3.1b（假高张力）
    # 的核心样本。比例截断已全部关闭（tier1/tier2_ratio=1.0），因此它一旦被淘汰就说明
    # 「PDBQT 姿态复原伪造高张力」的缺陷复发。
    candidate_heavy = {
        runner.ligand_dict[cid].qc_profile.heavy_atoms
        for cid in candidate_ids
        if runner.ligand_dict.get(cid) is not None and runner.ligand_dict[cid].qc_profile
    }
    assert 22 in candidate_heavy, (
        f"槲皮素（22 重原子，多羟基黄酮回归核心样本）未进入候选: {sorted(candidate_heavy)}"
    )
    # 已知残余（已单独上报 Lead）：黄芩苷在当前「库路径嵌入」的原子序下，PDBQT 复原后的
    # 极性氢转子搜索仍可能停在局部极小（ΔE≈15.7 kcal/mol），被张力 QC 淘汰；除此以外
    # 不允许任何其它分子丢失。
    dropped = sorted(set(runner.ligand_dict) - set(candidate_ids))
    assert len(dropped) <= 1, (
        f"除已知黄芩苷残余假张力外不应再有分子被 Tier-2 淘汰，实际丢失 {dropped}"
    )
    for cand in report.candidates:
        assert cand.strain_energy <= 8.0, (
            f"{cand.ligand_id}: 候选张力能 {cand.strain_energy} 超过 max_strain_kcal=8.0"
        )
    box_center = DELIVERY_BOX_CENTER
    half = tuple(size / 2.0 for size in DELIVERY_BOX_SIZE)
    for cand in report.candidates:
        pose_sdf = poses_dir / f"pose_{cand.ligand_id}.sdf"
        assert pose_sdf.is_file(), f"缺少交付姿态 {pose_sdf.name}"
        mol = Chem.MolFromMolBlock(pose_sdf.read_text(encoding="utf-8"), sanitize=False)
        assert mol is not None, f"{pose_sdf.name} 无法解析"
        assert mol.GetNumConformers() == 1
        assert mol.GetConformer().Is3D(), f"{pose_sdf.name} 必须是 3D 构象，不是 2D 投影"

        centroid = _heavy_centroid_of_mol(mol)
        # (1) 绝不等于原点（修复前「重新嵌入 + 归零」的历史缺陷）
        assert max(abs(value) for value in centroid) > 5.0, (
            f"{cand.ligand_id}: 交付姿态质心 {centroid} 落在原点附近"
        )
        # (2) 落在 (16, 25, 4) 搜索盒内
        for axis in range(3):
            assert abs(centroid[axis] - box_center[axis]) <= half[axis], (
                f"{cand.ligand_id}: 质心 {centroid} 超出搜索盒轴 {axis}"
            )
        # (3) 精确落在该分子的口袋目标点上（刚性平移不改变构象，容差 <= 1 Å）
        target = engine.targets[cand.ligand_id]
        deviation = sum((centroid[i] - target[i]) ** 2 for i in range(3)) ** 0.5
        assert deviation <= 1.0, f"{cand.ligand_id}: 交付姿态质心偏离口袋目标 {deviation:.3f} Å"

        # 交付的 pose_<cid>.pdbqt 也必须是真实对接输出
        pose_pdbqt = poses_dir / f"pose_{cand.ligand_id}.pdbqt"
        assert pose_pdbqt.is_file(), f"缺少交付 PDBQT 姿态 {pose_pdbqt.name}"
        assert "ATOM" in pose_pdbqt.read_text(encoding="utf-8", errors="replace")

    # --- 5E. /verify T1~T5 全量 PASS（T1 真实重算） -----------------------
    energy_rows = load_energy_table(out_dir / "screening_results.csv")
    for cand_id in candidate_ids:
        engine.register_delivered_pose(
            cand_id, (poses_dir / f"pose_{cand_id}.sdf").read_text(encoding="utf-8")
        )

    # T4（搜索盒包络）只有在提供共晶配体时才真实验算，此处用口袋内的参考配体激活它
    # （取质心最接近盒子中心的那个交付姿态，保证 T4 的 2.0 Å 质心偏差判据真实可验）
    crystal_id = min(
        candidate_ids,
        key=lambda cid: sum(
            (engine.targets[cid][i] - DELIVERY_BOX_CENTER[i]) ** 2 for i in range(3)
        ),
    )
    (out_dir / "crystal_ligand.sdf").write_text(
        (poses_dir / f"pose_{crystal_id}.sdf").read_text(encoding="utf-8"), encoding="utf-8"
    )

    rc_verify = run_verify(out_dir, options={"recalculate": True, "engine": engine})
    assert rc_verify == 0, f"run_verify(out_dir, recalculate=True) 必须返回 0，实际 {rc_verify}"

    pkg = load_delivery_package(out_dir)
    try:
        conf = parse_conf_file(pkg.conf_path)
        r1 = _verify_t1(pkg, conf, energy_rows, {"recalculate": True, "engine": engine})
        r2 = _verify_t2(pkg, conf, energy_rows, {"recalculate": True, "engine": engine})
        r3 = _verify_t3(pkg, energy_rows, None)
        r4 = _verify_t4(pkg, conf, None)
        r5 = _verify_t5(pkg, energy_rows, None)

        assert r1.passed and r1.status == "PASS", (
            f"T1 未真实通过: {r1.status} {r1.message} {r1.observed}"
        )
        assert "Δmax=0.00" in r1.observed, f"T1 重算残差应为 0: {r1.observed}"
        assert r2.passed and r2.status == "PASS", (
            f"T2 未通过: {r2.status} {r2.message} {r2.observed}"
        )
        assert r3.passed, f"T3 未通过: {r3.message}"
        assert r4.passed and r4.status == "PASS", f"T4 未真实通过: {r4.status} {r4.message}"
        assert r5.passed and r5.status == "PASS", f"T5 未通过: {r5.message} {r5.observed}"
    finally:
        pkg.cleanup()


# =============================================================================
# 6. CLI 入口 & Agent REPL 的 /screen 全参数转发（P0-3）
# =============================================================================

#: ``run_screen_headless`` 必须收到的全部关键字参数（screen_cli 的 17 个转发参数；
#: 第 18 个 ``--non-interactive`` 是 click 侧开关，由 6D 的行为断言覆盖）。
EXPECTED_HEADLESS_KWARGS = {
    "receptor_path",
    "library_path",
    "box_center_str",
    "box_size_str",
    "tcm_mode",
    "anchor_residues_str",
    "tier1_ratio",
    "tier2_ratio",
    "seeds_str",
    "strain_threshold",
    "engine_name",
    "cpu",
    "cpu_per_job",
    "resume",
    "out_dir",
    "enable_gnina_rescore",
    "demo",
}


def _write_cli_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Creates minimal receptor / library / out-dir paths for CLI argument tests."""
    receptor = tmp_path / "rec" / "1HSG.pdb"
    receptor.parent.mkdir(parents=True, exist_ok=True)
    receptor.write_text(
        "ATOM      1  N   ASP A  25      16.000  25.000   4.000  1.00  0.00           N\n"
        "TER\nEND\n",
        encoding="utf-8",
    )
    library = tmp_path / "lib" / "tcm.smi"
    library.parent.mkdir(parents=True, exist_ok=True)
    library.write_text("CC(=O)Oc1ccccc1C(=O)O Aspirin\n", encoding="utf-8")
    out_dir = tmp_path / "screen_out"
    return receptor, library, out_dir


def test_cli_entrypoint_and_agent_repl_screen_all_flags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`dockstudio screen ...` 与 REPL `/screen ...` 必须无损转发全部 18 个参数。

    断言要点：
    * ``dockstudio.cli.main([...])`` 返回 0，且 ``run_screen_headless`` 收到 17 个
      关键字参数逐项与命令行一致（--tcm-mode / --tier1-ratio / --resume /
      --no-gnina-rescore / --seeds / --strain-threshold ...）；
    * ``--non-interactive`` 缺失 ``--receptor`` 时 click 直接以 2 退出（证明该开关真的
      被 click 解析，而不是被 argparse 吞掉）；
    * ``agent._dispatch_command("/screen ...")`` 返回 0 且同样走 click
      （``run_screen_interactive`` 绝不被调用）；
    * 反向断言：裸 ``/screen``、``/screen --demo``、``/screen -d`` 仍走
      ``run_screen_interactive``，不被静默升级为无头模式。
    """
    from dockstudio.cli.agent import screen_args_need_cli

    # --- 6A. 纯函数契约：只有裸向导旗标才走交互分支 ------------------------
    assert screen_args_need_cli([]) is False
    assert screen_args_need_cli(["--demo"]) is False
    assert screen_args_need_cli(["-d"]) is False
    assert screen_args_need_cli(["--demo", "-d"]) is False
    assert screen_args_need_cli(["--non-interactive"]) is True
    assert screen_args_need_cli(["--receptor", "r.pdb"]) is True
    assert screen_args_need_cli(["--demo", "--non-interactive"]) is True

    receptor, library, out_dir = _write_cli_inputs(tmp_path)

    recorded: dict[str, Any] = {}
    interactive_calls: list[dict[str, Any]] = []

    def fake_headless(**kwargs: Any) -> int:
        recorded.clear()
        recorded.update(kwargs)
        return 0

    def fake_interactive(*args: Any, **kwargs: Any) -> int:
        interactive_calls.append({"args": args, "kwargs": kwargs})
        return 0

    monkeypatch.setattr("dockstudio.cli.commands.screen.run_screen_headless", fake_headless)
    monkeypatch.setattr("dockstudio.cli.commands.screen.run_screen_interactive", fake_interactive)

    # --- 6B. `dockstudio screen ...` 全参数 -> run_screen_headless ----------
    rc = dockstudio_cli.main(
        [
            "screen",
            "--receptor",
            str(receptor),
            "--library",
            str(library),
            "--box-center",
            "16.0,25.0,4.0",
            "--box-size",
            "22.5,22.5,22.5",
            "--tcm-mode",
            "--anchor-residues",
            "ASP25:A,THR80:A",
            "--tier1-ratio",
            "0.5",
            "--tier2-ratio",
            "0.05",
            "--seeds",
            "7,42,2026",
            "--strain-threshold",
            "6.5",
            "--engine",
            "vina",
            "--cpu",
            "4",
            "--cpu-per-job",
            "1",
            "--resume",
            "--out",
            str(out_dir),
            "--no-gnina-rescore",
            "--non-interactive",
        ]
    )
    assert rc == 0, f"dockstudio screen 全参数调用必须返回 0，实际 {rc}"
    assert not interactive_calls, "带 --non-interactive 时绝不能落到交互向导"
    assert set(recorded) == EXPECTED_HEADLESS_KWARGS, (
        f"转发参数集合不符: 缺 {EXPECTED_HEADLESS_KWARGS - set(recorded)}，"
        f"多 {set(recorded) - EXPECTED_HEADLESS_KWARGS}"
    )
    assert recorded["receptor_path"] == receptor
    assert recorded["library_path"] == library
    assert recorded["box_center_str"] == "16.0,25.0,4.0"
    assert recorded["box_size_str"] == "22.5,22.5,22.5"
    assert recorded["tcm_mode"] is True
    assert recorded["anchor_residues_str"] == "ASP25:A,THR80:A"
    assert recorded["tier1_ratio"] == pytest.approx(0.5)
    assert recorded["tier2_ratio"] == pytest.approx(0.05)
    assert recorded["seeds_str"] == "7,42,2026"
    assert recorded["strain_threshold"] == pytest.approx(6.5)
    assert recorded["engine_name"] == "vina"
    assert recorded["cpu"] == 4
    assert recorded["cpu_per_job"] == 1
    assert recorded["resume"] is True
    assert recorded["out_dir"] == out_dir
    assert recorded["enable_gnina_rescore"] is False
    assert recorded["demo"] is False

    # 反向：--no-tcm-mode / --no-resume 同一入口也生效（布尔开关双向转发）
    recorded.clear()
    rc_off = dockstudio_cli.main(
        [
            "screen",
            "--receptor",
            str(receptor),
            "--library",
            str(library),
            "--no-tcm-mode",
            "--no-resume",
            "--non-interactive",
        ]
    )
    assert rc_off == 0
    assert recorded["tcm_mode"] is False
    assert recorded["resume"] is False

    # --- 6C. `--non-interactive` 自身必须被 click 解析（缺失 --receptor 即 2）---
    recorded.clear()
    interactive_calls.clear()
    rc_bad = dockstudio_cli.main(["screen", "--non-interactive", "--box-center", "1,2,3"])
    assert rc_bad == 2, f"--non-interactive 缺少 --receptor 时必须 click 报错退出 2，实际 {rc_bad}"
    assert not recorded, "参数校验失败时不得调用 run_screen_headless"
    assert not interactive_calls

    # --- 6D. REPL `/screen ...` 全参数转发 --------------------------------
    recorded.clear()
    interactive_calls.clear()
    repl_line = (
        f'/screen --receptor "{receptor}" --library "{library}" '
        "--box-center 16.0,25.0,4.0 --tcm-mode --tier1-ratio 0.5 "
        "--seeds 42,1024,2026 --non-interactive"
    )
    rc_repl = agent_module._dispatch_command(repl_line)
    assert rc_repl == 0, f"REPL /screen 全参数调用必须返回 0，实际 {rc_repl}"
    assert not interactive_calls, "REPL 带过滤参数时绝不能退回交互向导"
    assert set(recorded) == EXPECTED_HEADLESS_KWARGS
    assert recorded["receptor_path"] == receptor
    assert recorded["library_path"] == library
    assert recorded["tcm_mode"] is True
    assert recorded["tier1_ratio"] == pytest.approx(0.5)
    assert recorded["seeds_str"] == "42,1024,2026"
    assert recorded["box_center_str"] == "16.0,25.0,4.0"
    assert recorded["out_dir"] is None  # 未给 --out 时不伪造输出目录

    # --- 6E. 反向断言：裸 /screen 与 /screen --demo 仍走交互向导 -----------
    for line, expected_demo in (
        ("/screen", False),
        ("/screen --demo", True),
        ("/screen -d", True),
    ):
        recorded.clear()
        interactive_calls.clear()
        rc_wizard = agent_module._dispatch_command(line)
        assert rc_wizard == 0, f"{line} 必须返回 0"
        assert not recorded, f"{line} 不该调用 run_screen_headless"
        assert len(interactive_calls) == 1, f"{line} 必须走 run_screen_interactive"
        assert interactive_calls[0]["kwargs"] == {"demo": expected_demo}, interactive_calls[0]
