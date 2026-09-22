#!/usr/bin/env python3
"""Run a dependency-free Morris screening for the server TECO-CNP model."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


TREATMENTS = (("P0", 0, 0.0), ("P25", 1, 25.0),
              ("P50", 1, 50.0), ("P100", 1, 100.0))


@dataclass(frozen=True)
class Parameter:
    parameter_id: str
    environment_variable: str
    parameter_name_cn: str
    primary_process: str
    source_in_model: str
    lower: float
    upper: float
    scale: str
    baseline_note: str

    def physical_value(self, unit_value: float) -> float:
        if self.scale == "log":
            return math.exp(math.log(self.lower) + unit_value *
                            (math.log(self.upper) - math.log(self.lower)))
        return self.lower + unit_value * (self.upper - self.lower)


def read_parameters(path: Path) -> List[Parameter]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty candidate table: {path}")
    result = []
    for row in rows:
        result.append(Parameter(
            parameter_id=row["parameter_id"].strip(),
            environment_variable=row["environment_variable"].strip(),
            parameter_name_cn=row["parameter_name_cn"].strip(),
            primary_process=row["primary_process"].strip(),
            source_in_model=row["source_in_model"].strip(),
            lower=float(row["lower_multiplier"]),
            upper=float(row["upper_multiplier"]),
            scale=row["scale"].strip().lower(),
            baseline_note=row["baseline_note"].strip(),
        ))
    return result


def write_csv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def morris_design(parameters: Sequence[Parameter], trajectories: int,
                  levels: int, seed: int) -> Tuple[List[dict], float]:
    if levels < 4 or levels % 2:
        raise ValueError("--levels must be an even integer >= 4")
    rng = random.Random(seed)
    step = 1.0 / (levels - 1)
    delta = levels / (2.0 * (levels - 1))
    grid = [i * step for i in range(levels)]
    low = [x for x in grid if x <= 1.0 - delta + 1e-12]
    high = [x for x in grid if x >= delta - 1e-12]
    design: List[dict] = []
    sample_number = 0

    for trajectory in range(1, trajectories + 1):
        directions = [rng.choice((-1, 1)) for _ in parameters]
        current = [rng.choice(low if d > 0 else high) for d in directions]
        order = list(range(len(parameters)))
        rng.shuffle(order)

        design.append({
            "sample_id": f"t{trajectory:03d}_s000",
            "sample_number": sample_number,
            "trajectory": trajectory,
            "step": 0,
            "changed_parameter": "",
            "unit_values": current.copy(),
        })
        sample_number += 1

        for step_number, index in enumerate(order, start=1):
            current = current.copy()
            current[index] += directions[index] * delta
            current[index] = min(1.0, max(0.0, current[index]))
            design.append({
                "sample_id": f"t{trajectory:03d}_s{step_number:03d}",
                "sample_number": sample_number,
                "trajectory": trajectory,
                "step": step_number,
                "changed_parameter": parameters[index].parameter_id,
                "unit_values": current.copy(),
            })
            sample_number += 1
    return design, delta


def parse_summary(path: Path) -> Dict[str, float]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected one summary row in {path}, found {len(rows)}")
    result = {key: float(value) for key, value in rows[0].items()}
    if not all(math.isfinite(value) for value in result.values()):
        raise RuntimeError(f"non-finite value in {path}")
    return result


def run_one(executable: Path, main_dir: Path, cache_dir: Path,
            parameters: Sequence[Parameter], sample_id: str,
            unit_values: Sequence[float], treatment: Tuple[str, int, float],
            start_year: int, end_year: int, timeout: int,
            force: bool = False) -> Tuple[str, str, Dict[str, float], float]:
    label, addition_flag, rate = treatment
    summary_path = (cache_dir / f"{sample_id}_{label}.csv").resolve()
    if summary_path.exists() and not force:
        return sample_id, label, parse_summary(summary_path), 0.0

    env = os.environ.copy()
    for parameter in parameters:
        env.pop(parameter.environment_variable, None)
    for parameter, unit_value in zip(parameters, unit_values):
        env[parameter.environment_variable] = f"{parameter.physical_value(unit_value):.12g}"
    env["TECO_SENSITIVITY_MODE"] = "1"
    env["TECO_SUMMARY_FILE"] = str(summary_path)

    command = [str(executable), str(start_year), str(end_year), "3", "0", "0", "1", "0",
               str(addition_flag), f"{rate:g}"]
    started = time.time()
    completed = subprocess.run(
        command, cwd=str(main_dir), env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=timeout, errors="replace",
    )
    elapsed = time.time() - started
    if completed.returncode != 0 or not summary_path.exists():
        failure_dir = cache_dir.parent / "failed_logs"
        failure_dir.mkdir(parents=True, exist_ok=True)
        log_path = failure_dir / f"{sample_id}_{label}.log"
        log_path.write_text(completed.stdout or "", encoding="utf-8")
        raise RuntimeError(
            f"TECO failed for {sample_id}/{label}, return code={completed.returncode}; "
            f"see {log_path}"
        )
    return sample_id, label, parse_summary(summary_path), elapsed


def case_values(treatment_metrics: Mapping[str, Mapping[str, float]]) -> Dict[str, float]:
    result: Dict[str, float] = {}
    p0 = treatment_metrics["P0"]
    for label, _, _ in TREATMENTS:
        for metric, value in treatment_metrics[label].items():
            result[f"{label}|{metric}"] = value
            if label != "P0":
                result[f"d{label}|{metric}"] = value - p0[metric]
    return result


def analyse_morris(design: Sequence[dict], delta: float,
                   metrics_by_sample: Mapping[str, Mapping[str, float]],
                   parameters: Sequence[Parameter], metric_groups: Mapping[str, Sequence[str]],
                   top_per_process: int, result_dir: Path) -> None:
    parameter_lookup = {item.parameter_id: item for item in parameters}
    effects: Dict[Tuple[str, str], List[float]] = {}
    by_trajectory: Dict[int, List[dict]] = {}
    for row in design:
        by_trajectory.setdefault(int(row["trajectory"]), []).append(row)

    for rows in by_trajectory.values():
        rows.sort(key=lambda item: int(item["step"]))
        for previous, current in zip(rows, rows[1:]):
            parameter_id = current["changed_parameter"]
            parameter_index = next(i for i, p in enumerate(parameters)
                                   if p.parameter_id == parameter_id)
            dx = current["unit_values"][parameter_index] - previous["unit_values"][parameter_index]
            if abs(dx) < 1e-12:
                raise RuntimeError(f"zero Morris step for {parameter_id}")
            before = metrics_by_sample[previous["sample_id"]]
            after = metrics_by_sample[current["sample_id"]]
            for output_name in before:
                effects.setdefault((parameter_id, output_name), []).append(
                    (after[output_name] - before[output_name]) / dx
                )

    ee_rows = []
    for (parameter_id, output_name), values in sorted(effects.items()):
        case, metric = output_name.split("|", 1)
        process = next((name for name, names in metric_groups.items() if metric in names), "")
        mean = statistics.fmean(values)
        mu_star = statistics.fmean(abs(value) for value in values)
        sigma = statistics.stdev(values) if len(values) > 1 else 0.0
        ee_rows.append({
            "parameter_id": parameter_id,
            "parameter_name_cn": parameter_lookup[parameter_id].parameter_name_cn,
            "primary_process": parameter_lookup[parameter_id].primary_process,
            "evaluated_process": process,
            "case": case,
            "metric": metric,
            "mu": mean,
            "mu_star": mu_star,
            "sigma": sigma,
            "n_effects": len(values),
        })
    write_csv(result_dir / "morris_elementary_effects.csv", ee_rows,
              ["parameter_id", "parameter_name_cn", "primary_process", "evaluated_process",
               "case", "metric", "mu", "mu_star", "sigma", "n_effects"])

    maxima: Dict[Tuple[str, str], float] = {}
    for row in ee_rows:
        if not row["evaluated_process"]:
            continue
        key = (str(row["case"]), str(row["metric"]))
        maxima[key] = max(maxima.get(key, 0.0), float(row["mu_star"]))

    process_scores = []
    for parameter in parameters:
        for process, metric_names in metric_groups.items():
            normalized = []
            for row in ee_rows:
                if row["parameter_id"] != parameter.parameter_id or row["metric"] not in metric_names:
                    continue
                maximum = maxima.get((str(row["case"]), str(row["metric"])), 0.0)
                if maximum > 0:
                    normalized.append(float(row["mu_star"]) / maximum)
            process_scores.append({
                "parameter_id": parameter.parameter_id,
                "parameter_name_cn": parameter.parameter_name_cn,
                "primary_process": parameter.primary_process,
                "evaluated_process": process,
                "process_score": statistics.fmean(normalized) if normalized else 0.0,
                "n_metric_cases": len(normalized),
            })

    for process in metric_groups:
        group = [row for row in process_scores if row["evaluated_process"] == process]
        group.sort(key=lambda row: float(row["process_score"]), reverse=True)
        for rank, row in enumerate(group, start=1):
            row["rank_within_process"] = rank
    write_csv(result_dir / "parameter_process_scores.csv", process_scores,
              ["parameter_id", "parameter_name_cn", "primary_process", "evaluated_process",
               "process_score", "n_metric_cases", "rank_within_process"])

    overall_rows = []
    for parameter in parameters:
        scores = [float(row["process_score"]) for row in process_scores
                  if row["parameter_id"] == parameter.parameter_id]
        own = next(float(row["process_score"]) for row in process_scores
                   if row["parameter_id"] == parameter.parameter_id and
                   row["evaluated_process"] == parameter.primary_process)
        overall_rows.append({
            "parameter_id": parameter.parameter_id,
            "parameter_name_cn": parameter.parameter_name_cn,
            "primary_process": parameter.primary_process,
            "own_process_score": own,
            "system_score": statistics.fmean(scores),
            "lower_multiplier": parameter.lower,
            "upper_multiplier": parameter.upper,
            "scale": parameter.scale,
        })
    overall_rows.sort(key=lambda row: float(row["system_score"]), reverse=True)
    for rank, row in enumerate(overall_rows, start=1):
        row["system_rank"] = rank
    write_csv(result_dir / "parameter_overall_ranking.csv", overall_rows,
              ["system_rank", "parameter_id", "parameter_name_cn", "primary_process",
               "own_process_score", "system_score", "lower_multiplier", "upper_multiplier", "scale"])

    selected = []
    for process in metric_groups:
        rows = [row for row in process_scores if row["evaluated_process"] == process]
        rows.sort(key=lambda row: float(row["process_score"]), reverse=True)
        for row in rows[:top_per_process]:
            selected.append(row)
    write_csv(result_dir / "selected_parameters_top_by_process.csv", selected,
              ["parameter_id", "parameter_name_cn", "primary_process", "evaluated_process",
               "process_score", "n_metric_cases", "rank_within_process"])


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="TECO-CNP Morris global sensitivity screening")
    parser.add_argument("--model-root", required=True, type=Path)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--candidates", type=Path,
                        default=script_dir / "teco_sensitivity_candidates.csv")
    parser.add_argument("--metric-groups", type=Path, default=script_dir / "metric_groups.json")
    parser.add_argument("--results", type=Path)
    parser.add_argument("--start-year", type=int, default=2021)
    parser.add_argument("--end-year", type=int, default=2024)
    parser.add_argument("--trajectories", type=int, default=20)
    parser.add_argument("--levels", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--top-per-process", type=int, default=3)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--design-only", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    args = parser.parse_args()

    model_root = args.model_root.expanduser().resolve()
    main_dir = model_root / "main"
    executable = (args.executable.expanduser().resolve() if args.executable else
                  main_dir / "teco_sensitivity.exe")
    result_dir = (args.results.expanduser().resolve() if args.results else
                  model_root / "output" / "sensitivity" / "TECO_C_P_coupling_Morris")
    result_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = result_dir / "run_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    parameters = read_parameters(args.candidates)
    metric_groups = json.loads(args.metric_groups.read_text(encoding="utf-8-sig"))
    design, delta = morris_design(parameters, args.trajectories, args.levels, args.seed)

    design_rows = []
    for row in design:
        out = {key: row[key] for key in
               ("sample_id", "sample_number", "trajectory", "step", "changed_parameter")}
        for parameter, unit_value in zip(parameters, row["unit_values"]):
            out[f"u_{parameter.parameter_id}"] = unit_value
            out[parameter.parameter_id] = parameter.physical_value(unit_value)
        design_rows.append(out)
    fields = list(design_rows[0].keys())
    write_csv(result_dir / "morris_design.csv", design_rows, fields)
    if args.design_only:
        print(f"Design written: {result_dir / 'morris_design.csv'}")
        return 0

    if not executable.is_file():
        raise FileNotFoundError(f"TECO executable not found: {executable}")

    # Baseline check: all multipliers equal one under all four P treatments.
    baseline_results = []
    baseline_unit = []
    for parameter in parameters:
        if parameter.scale == "log":
            baseline_unit.append((0.0 - math.log(parameter.lower)) /
                                 (math.log(parameter.upper) - math.log(parameter.lower)))
        else:
            baseline_unit.append((1.0 - parameter.lower) / (parameter.upper - parameter.lower))
    for treatment in TREATMENTS:
        _, label, metrics, elapsed = run_one(
            executable, main_dir, cache_dir, parameters, "baseline", baseline_unit,
            treatment, args.start_year, args.end_year, args.timeout, args.force)
        row = {"treatment": label, "elapsed_seconds": elapsed}
        row.update(metrics)
        baseline_results.append(row)
    write_csv(result_dir / "baseline_metrics.csv", baseline_results,
              list(baseline_results[0].keys()))
    if args.baseline_only:
        print(f"Baseline smoke test completed: {result_dir / 'baseline_metrics.csv'}")
        return 0

    tasks = []
    for row in design:
        for treatment in TREATMENTS:
            tasks.append((row, treatment))
    total = len(tasks)
    completed_count = 0
    all_results: Dict[str, Dict[str, Dict[str, float]]] = {}
    timings = []
    started_all = time.time()

    def execute(task):
        row, treatment = task
        return run_one(executable, main_dir, cache_dir, parameters,
                       row["sample_id"], row["unit_values"], treatment,
                       args.start_year, args.end_year, args.timeout, args.force)

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(execute, task): task for task in tasks}
        for future in as_completed(futures):
            sample_id, label, metrics, elapsed = future.result()
            all_results.setdefault(sample_id, {})[label] = metrics
            timings.append(elapsed)
            completed_count += 1
            if completed_count == 1 or completed_count % 20 == 0 or completed_count == total:
                elapsed_all = time.time() - started_all
                print(f"completed {completed_count}/{total} runs; elapsed {elapsed_all/60:.1f} min",
                      flush=True)

    long_rows = []
    metrics_by_sample: Dict[str, Dict[str, float]] = {}
    for design_row in design:
        sample_id = design_row["sample_id"]
        treatment_metrics = all_results[sample_id]
        metrics_by_sample[sample_id] = case_values(treatment_metrics)
        for label, _, _ in TREATMENTS:
            row = {"sample_id": sample_id, "treatment": label}
            row.update(treatment_metrics[label])
            long_rows.append(row)
    write_csv(result_dir / "sample_treatment_metrics.csv", long_rows,
              list(long_rows[0].keys()))

    analyse_morris(design, delta, metrics_by_sample, parameters,
                   metric_groups, args.top_per_process, result_dir)
    manifest = {
        "model_root": str(model_root),
        "executable": str(executable),
        "start_year": args.start_year,
        "end_year": args.end_year,
        "parameter_count": len(parameters),
        "trajectories": args.trajectories,
        "levels": args.levels,
        "delta": delta,
        "design_samples": len(design),
        "treatment_runs": total,
        "seed": args.seed,
        "workers": args.workers,
        "elapsed_seconds": time.time() - started_all,
        "mean_uncached_run_seconds": statistics.fmean([x for x in timings if x > 0])
        if any(x > 0 for x in timings) else 0.0,
    }
    (result_dir / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"TECO Morris screening completed: {result_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
