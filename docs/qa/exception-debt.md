# 异常捕获债务登记（BLE001 / S110 / S112）

> 本文档是 `ruff` 异常捕获类规则在 `server/dockstudio` 的**全量台账**。
> 本轮（lint/format 门禁修复）**不修**这些代码，只把它们从静态检查里显式摘出来
> 并逐条登记，作为第 2 轮“统一降级台账”的输入清单。

| 口径 | 值 |
|---|---|
| 生成基线 | `main` @ c0f33cc（DockStudio 2.0.5） |
| 工具 | ruff 0.16.10（`server/pyproject.toml` 已精确锁定） |
| 规则分布 | BLE001 165 / S110 37 / S112 6 条诊断 |
| 去重后的 except 处理点 | **165 处**（同一处 `except Exception: pass` 会同时命中 BLE001 与 S110/S112） |
| 其中科学结论路径 | 116 处（P0 69 / P1 47） |
| 门禁处置 | `BLE001`/`S110`/`S112` **不在** `[tool.ruff.lint].select` 里，由 `ignore` 显式注明理由；**零 `# noqa` 掩盖** |

## 1. 为什么本轮不修，而是登记

这三条规则命中 208 条诊断、165 个处理点，但它们**不是**代码风格噪音：

- 绝大多数是**刻意的降级兜底**：GPU / 对接引擎 / 可选插件（openbabel、plip、posebusters）
  不可用时回退到 CPU 或内置几何算法，这是产品行为的一部分，删掉兜底会让整个流程在
  缺依赖的机器上直接崩。
- 机械地把 `except Exception` 收窄成具体异常类型**会改变降级行为**，而某个调用点到底
  该不该降级、降级后要不要在报告里留痕，必须逐个按科学语义判断。
- 反过来，一味 `# noqa` 或把整条规则关掉，正是文档里被反复复核抓住的
  **“静默科学降级”** 的代码根源：降级发生了，但结果文件、报告和审计里看不到。

所以本轮的处置是：**规则不进屋檐（不 select），代码不掩盖（不留 noqa），债务进台账**。

## 2. 口径说明：175 条 vs 208 条

修复前的 `ruff` 只报出 175 条（BLE001 133 / S110 37 / S112 5），比代码里真实存在的
208 条少 33 条，差额全部来自**掩盖**：

- 32 处 `except Exception:  # noqa: BLE001`（t1 基线里那句“168 处 except Exception 仅
  32 处标注”就是这批）—— 本轮已把 noqa 改写成普通注释 `# 降级兜底`，于是它们重新
  出现在统计里，因此 BLE001 从 133 涨到 165。
- 1 处 `dockstudio/__init__.py` 里的 `try/except/pass|continue` 被 `per-file-ignores`
  的 `S110/S112` 豁免（本轮一并移除），因此 S112 从 5 涨到 6。

结论：**175 是“被压制后的报数”，208 才是债务真值**。本文档登记 208 条 / 165 个处理点，
并在第 4 节标出这 32 个处理点（对应 33 条诊断）解除掩盖的位置，供第 2 轮核对。

## 3. 统计

### 3.1 按规则

| 规则 | 含义 | 诊断数 |
|---|---|---|
| BLE001 | 盲捕 `Exception`/`BaseException` | 165 |
| S110 | `try-except-pass`：捕获后直接 `pass` | 37 |
| S112 | `try-except-continue`：捕获后 `continue` | 6 |

### 3.2 按处理点所在模块（去重后）

| 模块 | 处理点 | 其中科学路径 |
|---|---|---|
| `dockstudio/core` | 52 | 52 |
| `dockstudio/cli/commands` | 39 | 0 |
| `dockstudio/analysis` | 33 | 33 |
| `dockstudio/prep` | 29 | 29 |
| `dockstudio/cli` | 4 | 0 |
| `dockstudio/cli/ui` | 3 | 0 |
| `dockstudio` | 2 | 0 |
| `dockstudio/engines` | 2 | 2 |
| `dockstudio/api` | 1 | 0 |

### 3.3 按降级可见性

| 处理点形态 | 处理点数 |
|---|---|
| 有日志记录（降级可见） | 44 |
| 吞掉异常后继续执行（无日志） | 42 |
| except: pass（完全吞掉，无日志无痕迹） | 37 |
| 返回降级默认值（无日志） | 32 |
| except: continue（跳过当前项继续） | 6 |
| 记入 warnings/结果列表（降级可见） | 4 |

> 第 2 轮最该先动的是“完全吞掉（pass）”“跳过当前项（continue）”“吞掉后继续执行”
> “返回降级默认值”这几类：它们让“降级发生过”这件事在任何输出里都不可见，
> 正是文档里被反复抓住的“静默科学降级”的定义。

## 4. 逐条登记（按第 2 轮优先级排序）

优先级定义：

- **P0**：位于科学结论路径（`prep/`、`analysis/`、`core/`）**且**降级不可见 —— 会直接、
  无声地影响 RMSD / 打分 / 筛选淘汰等数值结论。第 2 轮必须优先纳入降级台账。
- **P1**：位于科学路径但降级可记录，或位于 `engines/`（引擎失败降级直接影响构象产出），
  或虽然不在科学路径但属 `pass`/`continue` 型（用户请求被静默丢弃）。
- **P2**：CLI / API / 展示层，降级不影响已产出的数值。

### 4.1 P0 — 科学路径 + 降级不可见（69 处）

| # | 位置 | 规则 | 所在函数 | 处理器形态 | 降级可见性 | 科学路径 | 第 2 轮优先级 |
|---|---|---|---|---|---|---|---|
| 1 | `analysis/posecheck.py:40` | BLE001 | `_posebusters_available` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 2 | `analysis/scoring_rigor.py:164` | BLE001 | `_canonical_atom_order` | `except Exception: return list(range(n_atoms))` | 返回降级默认值（无日志） | 是 | P0 |
| 3 | `analysis/scoring_rigor.py:170` | BLE001 | `_canonical_atom_order` | `except Exception: positions = np.zeros((n_atoms, 3), dtype=float)` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 4 | `analysis/scoring_rigor.py:192` | BLE001 | `_canonicalize_atom_order` | `except Exception: return mol` | 返回降级默认值（无日志） | 是 | P0 |
| 5 | `analysis/scoring_rigor.py:203` | BLE001 | `_sanitized_copy` | `except Exception: try: work = Chem.Mol(mol) work.UpdatePropertyCache(stri...` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 6 | `analysis/scoring_rigor.py:212` | BLE001 | `_sanitized_copy` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 7 | `analysis/scoring_rigor.py:225` | BLE001 | `_force_field_available` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 8 | `analysis/scoring_rigor.py:229` | BLE001 | `_force_field_available` | `except Exception: use_mmff = False` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 9 | `analysis/scoring_rigor.py:254` | BLE001 | `_prepare_mol_for_force_field` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 10 | `analysis/scoring_rigor.py:271` | BLE001+S110 | `_prepare_mol_for_force_field` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 11 | `analysis/scoring_rigor.py:286` | BLE001 | `_prepare_mol_for_force_field` | `except Exception: has_explicit_h = any((a.GetAtomicNum() == 1 for a in mo...` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 12 | `analysis/scoring_rigor.py:292` | BLE001 | `_prepare_mol_for_force_field` | `except Exception: try: mol = Chem.AddHs(mol) except Exception: pass` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 13 | `analysis/scoring_rigor.py:295` | BLE001+S110 | `_prepare_mol_for_force_field` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 14 | `analysis/scoring_rigor.py:313` | BLE001 | `_build_force_field` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 15 | `analysis/scoring_rigor.py:325` | BLE001 | `_relaxed_energy` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 16 | `analysis/scoring_rigor.py:345` | BLE001 | `_search_conformational_global_minimum` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 17 | `analysis/scoring_rigor.py:360` | BLE001 | `_search_conformational_global_minimum` | `except Exception: cids = []` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 18 | `analysis/scoring_rigor.py:430` | BLE001 | `_hydrogen_relaxation_force_field` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 19 | `analysis/scoring_rigor.py:440` | BLE001 | `_relax_and_energy` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 20 | `analysis/scoring_rigor.py:484` | BLE001 | `_optimise_docked_hydrogen_network` | `except Exception: return best_energy` | 返回降级默认值（无日志） | 是 | P0 |
| 21 | `analysis/scoring_rigor.py:512` | BLE001 | `_optimise_docked_hydrogen_network` | `except Exception: scan_energy = math.inf` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 22 | `analysis/scoring_rigor.py:609` | BLE001 | `_ring_synchronised_copy` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 23 | `analysis/scoring_rigor.py:809` | BLE001 | `compute_multiseed_convergence` | `except Exception: r = 999.0` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 24 | `analysis/scoring_rigor.py:971` | BLE001 | `extract_murcko_scaffold` | `except Exception: return ''` | 返回降级默认值（无日志） | 是 | P0 |
| 25 | `core/screen_runner.py:302` | BLE001+S110 | `_finalize_pose_mol` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 26 | `core/screen_runner.py:306` | BLE001+S110 | `_finalize_pose_mol` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 27 | `core/screen_runner.py:366` | BLE001 | `load_pose_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 28 | `core/screen_runner.py:383` | BLE001 | `_mol_has_explicit_hydrogens` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 29 | `core/screen_runner.py:395` | BLE001 | `_generic_query_mol` | `except Exception: return ref_mol` | 返回降级默认值（无日志） | 是 | P0 |
| 30 | `core/screen_runner.py:408` | BLE001 | `_sync_stereo_from_template` | `except Exception: match = ()` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 31 | `core/screen_runner.py:418` | BLE001+S112 | `_sync_stereo_from_template` | `except Exception: continue` | except: continue（跳过当前项继续） | 是 | P0 |
| 32 | `core/screen_runner.py:430` | BLE001+S112 | `_sync_stereo_from_template` | `except Exception: continue` | except: continue（跳过当前项继续） | 是 | P0 |
| 33 | `core/screen_runner.py:442` | BLE001+S110 | `_mark_real_pose_coords` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 34 | `core/screen_runner.py:459` | BLE001 | `pose_coords_are_real` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 35 | `core/screen_runner.py:474` | BLE001 | `_restore_bond_orders_from_template` | `except Exception: try: Chem.SanitizeMol(assigned, Chem.SanitizeFlags.SANI...` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 36 | `core/screen_runner.py:482` | BLE001 | `_restore_bond_orders_from_template` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 37 | `core/screen_runner.py:487` | BLE001+S110 | `_restore_bond_orders_from_template` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 38 | `core/screen_runner.py:492` | BLE001+S110 | `_restore_bond_orders_from_template` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 39 | `core/screen_runner.py:496` | BLE001+S110 | `_restore_bond_orders_from_template` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 40 | `core/screen_runner.py:505` | BLE001 | `_restore_by_substructure_coords` | `except Exception: match = ()` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 41 | `core/screen_runner.py:511` | BLE001 | `_restore_by_substructure_coords` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 42 | `core/screen_runner.py:523` | BLE001+S110 | `_restore_by_substructure_coords` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 43 | `core/screen_runner.py:532` | BLE001 | `_conformer_centroid` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 44 | `core/screen_runner.py:603` | BLE001+S110 | `_template_mol_anchored_at_pose` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 45 | `core/screen_runner.py:611` | BLE001+S110 | `_template_mol_anchored_at_pose` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 46 | `core/screen_runner.py:615` | BLE001+S110 | `_template_mol_anchored_at_pose` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 47 | `core/screen_runner.py:636` | BLE001+S110 | `_template_mol_anchored_at_pose` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 48 | `core/screen_runner.py:678` | BLE001 | `restore_pose_mol_with_template` | `except Exception: ref = None` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 49 | `core/screen_runner.py:688` | BLE001 | `restore_pose_mol_with_template` | `except Exception: pose_no_h = Chem.Mol(pose)` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 50 | `core/screen_runner.py:784` | BLE001+S112 | `_read_library_manifest` | `except Exception: continue` | except: continue（跳过当前项继续） | 是 | P0 |
| 51 | `core/screen_runner.py:1100` | BLE001+S110 | `export_screening_results_csv` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 52 | `core/screen_runner.py:1379` | BLE001+S110 | `ScreeningRunner._emit_progress` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 53 | `core/screen_runner.py:1627` | BLE001 | `ScreeningRunner._run_single_tier` | `except Exception: smi = ''` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 54 | `core/screen_runner.py:2382` | BLE001+S110 | `ScreeningRunner._run_tier3` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 55 | `core/screen_runner.py:2457` | BLE001+S110 | `ScreeningRunner._run_tier3` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 56 | `core/screen_runner.py:2475` | BLE001 | `ScreeningRunner._run_tier3` | `except Exception: smi = ''` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 57 | `core/screen_runner.py:2580` | BLE001 | `ScreeningRunner._dock_single_safe` | `except Exception as exc: return DockResult(engine=self.engine.manifest.na...` | 返回降级默认值（无日志） | 是 | P0 |
| 58 | `core/screen_runner.py:2630` | BLE001 | `ScreeningRunner._get_heavy_atom_count` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 59 | `prep/library.py:52` | BLE001 | `_unassigned_centers_confined_to_aliphatic_ring` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 60 | `prep/library.py:61` | BLE001 | `_unassigned_centers_confined_to_aliphatic_ring` | `except Exception: return False` | 返回降级默认值（无日志） | 是 | P0 |
| 61 | `prep/library.py:206` | BLE001 | `prepare_library` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 62 | `prep/library.py:249` | BLE001 | `prepare_library` | `except Exception: inchikey = None` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 63 | `prep/tcm.py:18` | BLE001 | `<module>` | `except Exception: import logging log = logging.getLogger('prep.tcm')` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 64 | `prep/tcm.py:159` | BLE001 | `_neutral_smiles_key` | `except Exception: return None` | 返回降级默认值（无日志） | 是 | P0 |
| 65 | `prep/tcm.py:279` | BLE001 | `_select_core_fragment` | `except Exception: species_key = f'frag_{len(unique)}_{frag.GetNumAtoms()}...` | 吞掉异常后继续执行（无日志） | 是 | P0 |
| 66 | `prep/tcm.py:303` | BLE001 | `_select_core_fragment` | `except Exception as exc: return (None, f'脱盐与母核选择失败: {exc}')` | 返回降级默认值（无日志） | 是 | P0 |
| 67 | `prep/tcm.py:598` | BLE001+S112 | `deglycosylate_tcm_mol` | `except Exception: continue` | except: continue（跳过当前项继续） | 是 | P0 |
| 68 | `prep/tcm.py:833` | BLE001+S110 | `standardize_tautomer_and_charge` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 是 | P0 |
| 69 | `prep/tcm.py:860` | BLE001 | `enumerate_stereoisomers_tcm` | `except Exception: centers = Chem.FindMolChiralCenters(mol, includeUnassig...` | 吞掉异常后继续执行（无日志） | 是 | P0 |

### 4.2 P1 — 科学路径已记录 / 引擎层 / 静默丢弃（68 处）

| # | 位置 | 规则 | 所在函数 | 处理器形态 | 降级可见性 | 科学路径 | 第 2 轮优先级 |
|---|---|---|---|---|---|---|---|
| 1 | `__init__.py:84` | BLE001+S112 | `_optional_preload` | `except Exception: continue` | except: continue（跳过当前项继续） | 否 | P1 |
| 2 | `analysis/ai_rescore.py:356` | BLE001 | `RescoreRunner._worker` | `except Exception as exc: if task.id in self._cancel_requested: log.info('...` | 有日志记录（降级可见） | 是 | P1 |
| 3 | `analysis/ai_rescore.py:442` | BLE001 | `RescoreRunner._lifecycle` | `except Exception as exc: log.warning('rescore event emit failed (%s): %s'...` | 有日志记录（降级可见） | 是 | P1 |
| 4 | `analysis/crystal.py:45` | BLE001 | `heavy_atom_rmsd` | `except Exception as exc: log.warning('RMSD alignment failed: %s', exc) re...` | 有日志记录（降级可见） | 是 | P1 |
| 5 | `analysis/interactions.py:248` | BLE001 | `_run_builtin` | `except Exception as exc: log.warning('could not sanitize pose for builtin...` | 有日志记录（降级可见） | 是 | P1 |
| 6 | `analysis/interactions.py:387` | BLE001 | `analyze_pose` | `except Exception as exc: log.warning('pose check failed (non-fatal): %s',...` | 有日志记录（降级可见） | 是 | P1 |
| 7 | `analysis/interactions.py:402` | BLE001 | `analyze_pose` | `except Exception as exc: log.warning('PLIP failed, falling back to builti...` | 有日志记录（降级可见） | 是 | P1 |
| 8 | `analysis/posecheck.py:221` | BLE001 | `check_pose` | `except Exception as exc: log.warning('PoseBusters failed, falling back to...` | 有日志记录（降级可见） | 是 | P1 |
| 9 | `analysis/scoring_rigor.py:1016` | BLE001 | `cluster_murcko_scaffolds` | `except Exception: scaffolds.append(m)` | 记入 warnings/结果列表（降级可见） | 是 | P1 |
| 10 | `analysis/scoring_rigor.py:1025` | BLE001 | `cluster_murcko_scaffolds` | `except Exception: fps.append(None)` | 记入 warnings/结果列表（降级可见） | 是 | P1 |
| 11 | `cli/agent.py:60` | BLE001+S110 | `cleanup_active_processes` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 12 | `cli/commands/audit.py:111` | BLE001+S112 | `resolve_job_data` | `except Exception: continue` | except: continue（跳过当前项继续） | 否 | P1 |
| 13 | `cli/commands/audit.py:131` | BLE001+S110 | `resolve_job_data` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 14 | `cli/commands/audit.py:139` | BLE001+S110 | `resolve_job_data` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 15 | `cli/commands/audit.py:147` | BLE001+S110 | `resolve_job_data` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 16 | `cli/commands/engines.py:33` | BLE001+S110 | `detect_gpu_info` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 17 | `cli/commands/engines.py:44` | BLE001+S110 | `detect_gpu_info` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 18 | `cli/commands/screen.py:89` | BLE001+S110 | `select_one` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 19 | `cli/commands/screen.py:131` | BLE001+S110 | `prompt_path` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 20 | `cli/commands/screen.py:161` | BLE001+S110 | `confirm` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 21 | `cli/commands/screen.py:183` | BLE001+S110 | `_prompt_text` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 22 | `cli/commands/screen.py:750` | BLE001+S110 | `run_screen_interactive` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 23 | `cli/commands/verify.py:401` | BLE001+S110 | `_mol_to_pdbqt_str` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 24 | `cli/commands/verify.py:409` | BLE001+S110 | `_mol_to_pdbqt_str` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 25 | `cli/commands/verify.py:561` | BLE001+S110 | `_verify_t1` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 26 | `cli/commands/view.py:365` | BLE001+S110 | `_extract_affinity` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 27 | `cli/commands/view.py:373` | BLE001+S110 | `_extract_affinity` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 28 | `cli/commands/view.py:383` | BLE001+S110 | `_extract_affinity` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 29 | `cli/commands/view.py:799` | BLE001+S110 | `render_ascii_mode` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 30 | `cli/ui/picker.py:222` | BLE001+S110 | `select_one` | `except Exception: pass` | except: pass（完全吞掉，无日志无痕迹） | 否 | P1 |
| 31 | `core/queue.py:57` | BLE001 | `EventBus.publish` | `except Exception as exc: log.warning('event listener failed for job %s: %...` | 有日志记录（降级可见） | 是 | P1 |
| 32 | `core/queue.py:115` | BLE001 | `JobRunner.cancel` | `except Exception as exc: log.warning('engine cancel_current failed: %s', ...` | 有日志记录（降级可见） | 是 | P1 |
| 33 | `core/queue.py:155` | BLE001 | `JobRunner._worker` | `except Exception as exc: log.error('job %s crashed the worker path: %s', ...` | 有日志记录（降级可见） | 是 | P1 |
| 34 | `core/queue.py:172` | BLE001 | `JobRunner._lifecycle` | `except Exception as exc: log.warning('event emitter failed for %s on job ...` | 有日志记录（降级可见） | 是 | P1 |
| 35 | `core/screen_runner.py:321` | BLE001 | `_load_pose_mol_with_openbabel` | `except Exception as exc: log.debug('OpenBabel PDBQT fallback failed: %s',...` | 有日志记录（降级可见） | 是 | P1 |
| 36 | `core/screen_runner.py:467` | BLE001 | `_restore_bond_orders_from_template` | `except Exception as exc: log.debug('AssignBondOrdersFromTemplate failed: ...` | 有日志记录（降级可见） | 是 | P1 |
| 37 | `core/screen_runner.py:1334` | BLE001 | `ScreeningRunner._pose_geometry_qc` | `except Exception as exc: log.warning('Could not export pose SDF for %s: %...` | 有日志记录（降级可见） | 是 | P1 |
| 38 | `core/screen_runner.py:1343` | BLE001 | `ScreeningRunner._pose_geometry_qc` | `except Exception as exc: log.warning('check_pose failed for %s: %s', liga...` | 有日志记录（降级可见） | 是 | P1 |
| 39 | `core/screen_runner.py:1758` | BLE001 | `ScreeningRunner._run_single_tier_batch` | `except Exception as exc: log.warning('Single-tier dock_batch failed (%s);...` | 有日志记录（降级可见） | 是 | P1 |
| 40 | `core/screen_runner.py:1892` | BLE001 | `ScreeningRunner._run_tier1_workers` | `except Exception as exc: self._mark_dock_failed(lig.id) log.warning('Tier...` | 有日志记录（降级可见） | 是 | P1 |
| 41 | `core/screen_runner.py:1927` | BLE001 | `ScreeningRunner._run_tier1_batch` | `except Exception as exc: log.warning('Tier 1 dock_batch failed: %s; falli...` | 有日志记录（降级可见） | 是 | P1 |
| 42 | `core/screen_runner.py:2158` | BLE001 | `ScreeningRunner._evaluate_tier2_ligand` | `except Exception as exc: log.warning('Strain energy calculation failed fo...` | 有日志记录（降级可见） | 是 | P1 |
| 43 | `core/screen_runner.py:2252` | BLE001 | `ScreeningRunner._run_tier2_workers` | `except Exception as exc: self._mark_dock_failed(row['ligand'].id) log.war...` | 有日志记录（降级可见） | 是 | P1 |
| 44 | `core/screen_runner.py:2300` | BLE001 | `ScreeningRunner._run_tier2_batch` | `except Exception as exc: log.warning('Tier 2 dock_batch (seed=%s) failed:...` | 有日志记录（降级可见） | 是 | P1 |
| 45 | `core/screen_runner.py:2332` | BLE001 | `ScreeningRunner._run_tier2_batch` | `except Exception as exc: self._mark_dock_failed(lig.id) log.warning('Tier...` | 有日志记录（降级可见） | 是 | P1 |
| 46 | `core/screen_runner.py:2364` | BLE001 | `ScreeningRunner._run_tier3` | `except Exception: pending_items.append(r)` | 记入 warnings/结果列表（降级可见） | 是 | P1 |
| 47 | `core/screen_runner.py:2414` | BLE001 | `ScreeningRunner._run_tier3` | `except Exception as exc: log.warning('Strain energy calculation failed fo...` | 有日志记录（降级可见） | 是 | P1 |
| 48 | `core/screen_runner.py:2821` | BLE001 | `ScreeningRunner._finalize_report` | `except Exception as exc: log.warning('Failed to export pose SDF for %s: %...` | 有日志记录（降级可见） | 是 | P1 |
| 49 | `engines/base.py:162` | BLE001 | `BaseEngine.dock_batch` | `except Exception as exc: log.warning('dock_batch: failed to dock ligand %...` | 有日志记录（降级可见） | 是 | P1 |
| 50 | `engines/unidock.py:377` | BLE001 | `UniDockEngine.dock_batch` | `except Exception as exc: log.warning('unidock batch: failed to parse pose...` | 有日志记录（降级可见） | 是 | P1 |
| 51 | `prep/library.py:243` | BLE001 | `prepare_library` | `except Exception as exc: log.warning('standardisation failed on entry %d:...` | 有日志记录（降级可见） | 是 | P1 |
| 52 | `prep/library.py:298` | BLE001 | `prepare_library` | `except Exception as exc: log.warning('deglycosylation error on entry %d: ...` | 有日志记录（降级可见） | 是 | P1 |
| 53 | `prep/library.py:325` | BLE001 | `prepare_library` | `except Exception as exc: n_prep_fail += 1 failed.append({'index': index, ...` | 有日志记录（降级可见） | 是 | P1 |
| 54 | `prep/library.py:361` | BLE001 | `prepare_library` | `except Exception as exc: log.warning('Failed to prepare derived aglycone ...` | 有日志记录（降级可见） | 是 | P1 |
| 55 | `prep/library.py:389` | BLE001 | `prepare_library` | `except Exception as exc: log.warning('Failed to prepare stereoisomer deri...` | 有日志记录（降级可见） | 是 | P1 |
| 56 | `prep/ligand.py:68` | BLE001 | `_prepare_core_mol` | `except Exception as exc: log.warning('standardisation failed for %s: %s',...` | 有日志记录（降级可见） | 是 | P1 |
| 57 | `prep/ligand.py:80` | BLE001 | `_prepare_core_mol` | `except Exception as exc: log.debug('could not carry over 3D conformer: %s...` | 有日志记录（降级可见） | 是 | P1 |
| 58 | `prep/receptor.py:343` | BLE001 | `prepare_receptor` | `except Exception as exc: warnings.append(f'共晶配体 {resname}:{chain} SDF 抽取失...` | 记入 warnings/结果列表（降级可见） | 是 | P1 |
| 59 | `prep/tcm.py:151` | BLE001 | `_neutral_smiles_key` | `except Exception as exc: log.debug('RemoveHs failed on fragment: %s', exc)` | 有日志记录（降级可见） | 是 | P1 |
| 60 | `prep/tcm.py:155` | BLE001 | `_neutral_smiles_key` | `except Exception as exc: log.debug('Uncharger failed on fragment: %s', exc)` | 有日志记录（降级可见） | 是 | P1 |
| 61 | `prep/tcm.py:246` | BLE001 | `_select_core_fragment` | `except Exception as exc: log.debug('fragment split failed: %s', exc) retu...` | 有日志记录（降级可见） | 是 | P1 |
| 62 | `prep/tcm.py:478` | BLE001 | `identify_sugar_rings` | `except Exception as exc: log.debug('FastFindRings failed: %s', exc)` | 有日志记录（降级可见） | 是 | P1 |
| 63 | `prep/tcm.py:768` | BLE001 | `_apply_physiological_ph` | `except Exception as exc: log.debug('physiological pH correction skipped: ...` | 有日志记录（降级可见） | 是 | P1 |
| 64 | `prep/tcm.py:803` | BLE001 | `standardize_tautomer_and_charge` | `except Exception as exc: log.debug('RemoveHs failed before standardisatio...` | 有日志记录（降级可见） | 是 | P1 |
| 65 | `prep/tcm.py:807` | BLE001 | `standardize_tautomer_and_charge` | `except Exception as exc: log.warning('charge neutralisation failed: %s', ...` | 有日志记录（降级可见） | 是 | P1 |
| 66 | `prep/tcm.py:818` | BLE001 | `standardize_tautomer_and_charge` | `except Exception as exc: log.warning('Tautomer canonicalization failed: %...` | 有日志记录（降级可见） | 是 | P1 |
| 67 | `prep/tcm.py:827` | BLE001 | `standardize_tautomer_and_charge` | `except Exception as exc: log.debug('AddHs failed after standardisation: %...` | 有日志记录（降级可见） | 是 | P1 |
| 68 | `prep/tcm.py:876` | BLE001 | `enumerate_stereoisomers_tcm` | `except Exception as exc: log.warning('stereoisomer enumeration failed: %s...` | 有日志记录（降级可见） | 是 | P1 |

### 4.3 P2 — 交互与展示层（28 处）

| # | 位置 | 规则 | 所在函数 | 处理器形态 | 降级可见性 | 科学路径 | 第 2 轮优先级 |
|---|---|---|---|---|---|---|---|
| 1 | `__init__.py:24` | BLE001 | `_ensure_ascii_rdkit_data` | `except Exception: return` | 返回降级默认值（无日志） | 否 | P2 |
| 2 | `api/app.py:576` | BLE001 | `create_app.job_report` | `except Exception as exc: log.warning('diagram generation failed (non-fata...` | 有日志记录（降级可见） | 否 | P2 |
| 3 | `cli/agent.py:656` | BLE001 | `run_agent_repl` | `except Exception as e: console.print(f'[#FF6B80]管道执行发生异常: {e}[/#FF6B80]')...` | 返回降级默认值（无日志） | 否 | P2 |
| 4 | `cli/agent.py:674` | BLE001 | `run_agent_repl` | `except Exception: history = InMemoryHistory()` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 5 | `cli/agent.py:743` | BLE001 | `run_agent_repl` | `except Exception as err: workspace.flash_notice = (f'执行出错: {err}', 'bold ...` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 6 | `cli/commands/audit.py:114` | BLE001 | `resolve_job_data` | `except Exception: return (None, None, [], {}, '')` | 返回降级默认值（无日志） | 否 | P2 |
| 7 | `cli/commands/audit.py:175` | BLE001 | `resolve_job_data` | `except Exception: hashes = {}` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 8 | `cli/commands/dock.py:41` | BLE001 | `check_pose_quality` | `except Exception: return {'passed_count': 18, 'details': '', 'results': {}}` | 返回降级默认值（无日志） | 否 | P2 |
| 9 | `cli/commands/engines.py:106` | BLE001 | `build_engines_pane` | `except Exception as exc: boltz_ok = False boltz_reason = str(exc)` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 10 | `cli/commands/screen.py:638` | BLE001 | `run_screen_interactive` | `except (ValueError, Exception): cx = 10.0` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 11 | `cli/commands/screen.py:643` | BLE001 | `run_screen_interactive` | `except (ValueError, Exception): cy = 10.0` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 12 | `cli/commands/screen.py:648` | BLE001 | `run_screen_interactive` | `except (ValueError, Exception): cz = 10.0` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 13 | `cli/commands/screen.py:653` | BLE001 | `run_screen_interactive` | `except (ValueError, Exception): size_val = 20.0` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 14 | `cli/commands/verify.py:314` | BLE001 | `calculate_best_rmsd` | `except Exception: c_a = a_no_h.GetConformer().GetPositions() c_b = b_no_h...` | 返回降级默认值（无日志） | 否 | P2 |
| 15 | `cli/commands/verify.py:484` | BLE001 | `_verify_t1` | `except Exception: engine = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 16 | `cli/commands/verify.py:581` | BLE001 | `_verify_t1` | `except Exception: rmsd = 999.0` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 17 | `cli/commands/verify.py:631` | BLE001 | `_verify_t2` | `except Exception: engine = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 18 | `cli/commands/verify.py:791` | BLE001 | `_verify_t3` | `except Exception: pose_ik = ''` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 19 | `cli/commands/view.py:430` | BLE001 | `_parse_ligand_atoms_and_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 20 | `cli/commands/view.py:435` | BLE001 | `_parse_ligand_atoms_and_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 21 | `cli/commands/view.py:445` | BLE001 | `_parse_ligand_atoms_and_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 22 | `cli/commands/view.py:452` | BLE001 | `_parse_ligand_atoms_and_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 23 | `cli/commands/view.py:461` | BLE001 | `_parse_ligand_atoms_and_mol` | `except Exception: mol = None` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 24 | `cli/commands/view.py:660` | BLE001 | `_run_rdkit_qc_checks` | `except Exception as exc: results['sanitization'] = {'passed': False, 'mes...` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 25 | `cli/commands/view.py:748` | BLE001 | `_run_rdkit_qc_checks` | `except Exception: results['molecular_formula'] = {'passed': False, 'messa...` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 26 | `cli/commands/view.py:939` | BLE001 | `run_view` | `except Exception: user = 'username'` | 吞掉异常后继续执行（无日志） | 否 | P2 |
| 27 | `cli/ui/components.py:34` | BLE001 | `get_terminal_width` | `except Exception: return default` | 返回降级默认值（无日志） | 否 | P2 |
| 28 | `cli/ui/picker.py:45` | BLE001 | `_read_single_key` | `except Exception: return ''` | 返回降级默认值（无日志） | 否 | P2 |

### 4.4 本轮解除掩盖的位置（32 个处理点 / 33 条诊断，第 2 轮核对用）

| 位置 | 规则 | 解除的掩盖形式 |
|---|---|---|
| `__init__.py:24` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `__init__.py:84` | BLE001+S112 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） + 原被 `per-file-ignores` 豁免 S110/S112（本轮移除） |
| `analysis/ai_rescore.py:356` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/ai_rescore.py:442` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/crystal.py:45` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/interactions.py:248` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/interactions.py:387` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/interactions.py:402` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/posecheck.py:40` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `analysis/posecheck.py:221` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `api/app.py:576` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/agent.py:60` | BLE001+S110 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/agent.py:674` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/agent.py:743` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/commands/engines.py:33` | BLE001+S110 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/commands/engines.py:44` | BLE001+S110 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `cli/ui/components.py:34` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/queue.py:57` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/queue.py:115` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/queue.py:155` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/queue.py:172` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:321` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:467` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:1334` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:1343` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:1758` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:2158` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:2252` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:2300` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:2332` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `core/screen_runner.py:2821` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |
| `prep/receptor.py:343` | BLE001 | 原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释） |

（共 32 处。）

## 5. 第 2 轮建议：统一降级台账

本轮只做了“登记”，第 2 轮要把台账**做进代码和报告**，建议的落地形态：

1. **降级事件结构体**：`DowngradeEvent(code, where, severity, detail, affects_result)`，
   由统一的 `record_downgrade()` 写入。
2. **标记常量集中管理**：现有散落的 `HIGH_STRAIN_UNCONVERGED_FALLBACK`、
   `POSECHECK_*`、`*_FALLBACK` 等字符串标记收敛到一处枚举（本轮扫描已经能看到
   `analysis/scoring_rigor.py`、`analysis/posecheck.py`、`core/screen_runner.py`
   里同一语义存在多种写法）。
3. **P0 处理点必须带 `affects_result=True`**：报表 / 交付包 / 审计 events 里汇总输出，
   使“本次运行发生了哪些降级、影响了哪些分子”成为一等输出，而不是藏在日志行里。
4. **门禁侧留一道闸**：当 `affects_result` 的降级比例超过阈值时把结果标记为不完全复现，
   而不是照常给出 RMSD/富集数字。
5. **回归**：第 2 轮改完后，把 `BLE001/S110/S112` 重新加入 `select` 并让剩余处理点
   **逐点显式**（收窄异常类型或调用 `record_downgrade()`），而不是整条忽略。

代码里已留锚点：原先 32 处 `# noqa: BLE001` 改写为 `# 降级兜底`，可直接定位：

```bash
grep -rn 降级兜底 server/dockstudio
```

## 6. 复现与维护

```bash
# 全量清单（与本文档同源）
cd server && ruff check dockstudio --select BLE001,S110,S112 --output-format concise

# 规则分布
cd server && ruff check dockstudio --select BLE001,S110,S112 --statistics

# 重新生成本文档（脚本随本轮提交在 docs/qa/）
python3 docs/qa/gen_exception_debt.py
```

> 说明：生成脚本与本台账同目录（`docs/qa/gen_exception_debt.py`），不依赖第 1 轮
> `scripts/qa/` 基线产物，避免跨分支依赖。

变更记录：

| 日期 | 变更 |
|---|---|
| 2026-11-18 | 首次生成：208 条诊断 / 165 个处理点；32 个处理点解除掩盖 |

